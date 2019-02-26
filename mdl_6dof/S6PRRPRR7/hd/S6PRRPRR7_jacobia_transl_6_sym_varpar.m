% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:36
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (225->42), mult. (461->74), div. (0->0), fcn. (581->12), ass. (0->34)
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t23 = qJ(5) + qJ(6);
t21 = sin(t23);
t22 = cos(t23);
t26 = sin(qJ(5));
t34 = pkin(5) * t26 + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(4);
t39 = pkin(3) + r_i_i_C(3) + pkin(10) + pkin(9);
t48 = t34 * t27 + t39 * t30 + pkin(2);
t24 = sin(pkin(11));
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t40 = cos(pkin(11));
t41 = cos(pkin(6));
t36 = t41 * t40;
t11 = t24 * t28 - t31 * t36;
t12 = t24 * t31 + t28 * t36;
t25 = sin(pkin(6));
t37 = t25 * t40;
t7 = t12 * t27 + t30 * t37;
t47 = (-t11 * t21 + t7 * t22) * r_i_i_C(1) + (-t11 * t22 - t7 * t21) * r_i_i_C(2);
t38 = t24 * t41;
t13 = t40 * t28 + t31 * t38;
t14 = -t28 * t38 + t40 * t31;
t43 = t25 * t30;
t9 = t14 * t27 - t24 * t43;
t46 = (-t13 * t21 + t9 * t22) * r_i_i_C(1) + (-t13 * t22 - t9 * t21) * r_i_i_C(2);
t44 = t25 * t27;
t15 = t28 * t44 - t41 * t30;
t42 = t25 * t31;
t45 = (t15 * t22 + t21 * t42) * r_i_i_C(1) + (-t15 * t21 + t22 * t42) * r_i_i_C(2);
t29 = cos(qJ(5));
t35 = t29 * pkin(5) + t22 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(4) + pkin(8);
t1 = [0, -t13 * t48 + t35 * t14, -t39 * t9 + t34 * (t14 * t30 + t24 * t44) t9 (-t13 * t26 + t29 * t9) * pkin(5) + t46, t46; 0, -t11 * t48 + t35 * t12, -t39 * t7 + t34 * (t12 * t30 - t27 * t37) t7 (-t11 * t26 + t29 * t7) * pkin(5) + t47, t47; 1 (t35 * t28 + t48 * t31) * t25, -t39 * t15 + t34 * (t41 * t27 + t28 * t43) t15 (t15 * t29 + t26 * t42) * pkin(5) + t45, t45;];
Ja_transl  = t1;
