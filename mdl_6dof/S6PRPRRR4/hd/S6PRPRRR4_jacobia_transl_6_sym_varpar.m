% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (281->45), mult. (394->75), div. (0->0), fcn. (497->13), ass. (0->35)
t24 = pkin(12) + qJ(4);
t20 = sin(t24);
t21 = cos(t24);
t25 = qJ(5) + qJ(6);
t22 = sin(t25);
t23 = cos(t25);
t31 = cos(qJ(5));
t36 = t31 * pkin(5) + r_i_i_C(1) * t23 - r_i_i_C(2) * t22 + pkin(4);
t45 = r_i_i_C(3) + pkin(10) + pkin(9);
t49 = t45 * t20 + t36 * t21 + cos(pkin(12)) * pkin(3) + pkin(2);
t26 = sin(pkin(11));
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t40 = cos(pkin(11));
t41 = cos(pkin(6));
t37 = t41 * t40;
t13 = t26 * t30 - t32 * t37;
t14 = t26 * t32 + t30 * t37;
t27 = sin(pkin(6));
t38 = t27 * t40;
t8 = t14 * t21 - t20 * t38;
t48 = (t13 * t23 - t8 * t22) * r_i_i_C(1) + (-t13 * t22 - t8 * t23) * r_i_i_C(2);
t39 = t26 * t41;
t16 = -t30 * t39 + t40 * t32;
t44 = t26 * t27;
t10 = t16 * t21 + t20 * t44;
t15 = t40 * t30 + t32 * t39;
t47 = (-t10 * t22 + t15 * t23) * r_i_i_C(1) + (-t10 * t23 - t15 * t22) * r_i_i_C(2);
t43 = t27 * t30;
t12 = t41 * t20 + t21 * t43;
t42 = t27 * t32;
t46 = (-t12 * t22 - t23 * t42) * r_i_i_C(1) + (-t12 * t23 + t22 * t42) * r_i_i_C(2);
t29 = sin(qJ(5));
t35 = t29 * pkin(5) + t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(8) + qJ(3);
t1 = [0, -t15 * t49 + t35 * t16, t15, t45 * t10 + t36 * (-t16 * t20 + t21 * t44) (-t10 * t29 + t15 * t31) * pkin(5) + t47, t47; 0, -t13 * t49 + t35 * t14, t13, t45 * t8 + t36 * (-t14 * t20 - t21 * t38) (t13 * t31 - t29 * t8) * pkin(5) + t48, t48; 1 (t35 * t30 + t32 * t49) * t27, -t42, t45 * t12 + t36 * (-t20 * t43 + t41 * t21) (-t12 * t29 - t31 * t42) * pkin(5) + t46, t46;];
Ja_transl  = t1;
