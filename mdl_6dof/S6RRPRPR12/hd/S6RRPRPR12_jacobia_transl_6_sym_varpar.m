% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:22
% EndTime: 2019-02-26 21:44:22
% DurationCPUTime: 0.15s
% Computational Cost: add. (267->56), mult. (463->89), div. (0->0), fcn. (584->12), ass. (0->38)
t47 = pkin(10) + r_i_i_C(3);
t46 = pkin(2) + qJ(5) + pkin(9);
t23 = sin(pkin(6));
t27 = sin(qJ(2));
t45 = t23 * t27;
t28 = sin(qJ(1));
t44 = t23 * t28;
t31 = cos(qJ(2));
t43 = t23 * t31;
t32 = cos(qJ(1));
t42 = t23 * t32;
t41 = cos(pkin(6));
t30 = cos(qJ(4));
t40 = t23 * (t30 * pkin(4) + pkin(3) + pkin(8));
t26 = sin(qJ(4));
t39 = t26 * pkin(4) + qJ(3);
t38 = t28 * t41;
t37 = t32 * t41;
t25 = sin(qJ(6));
t29 = cos(qJ(6));
t36 = t29 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(5);
t13 = t28 * t27 - t31 * t37;
t22 = qJ(4) + pkin(11);
t20 = sin(t22);
t21 = cos(t22);
t7 = -t13 * t20 + t21 * t42;
t35 = t13 * t21 + t20 * t42;
t34 = t25 * r_i_i_C(1) + t29 * r_i_i_C(2) + t46;
t33 = t36 * t20 - t47 * t21 + t39;
t16 = -t27 * t38 + t32 * t31;
t15 = t32 * t27 + t31 * t38;
t14 = t27 * t37 + t28 * t31;
t12 = -t20 * t43 + t41 * t21;
t4 = t15 * t20 + t21 * t44;
t3 = -t15 * t21 + t20 * t44;
t2 = t16 * t25 + t4 * t29;
t1 = t16 * t29 - t4 * t25;
t5 = [-t28 * pkin(1) - t39 * t13 - t34 * t14 + t32 * t40 + t47 * t35 + t36 * t7, -t34 * t15 + t33 * t16, t15, t47 * t4 + (t15 * t30 - t26 * t44) * pkin(4) - t36 * t3, t16, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t32 * pkin(1) + t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * t15 + t46 * t16 + t28 * t40 + t47 * t3, -t34 * t13 + t33 * t14, t13, -t47 * t7 + (t13 * t30 + t26 * t42) * pkin(4) + t36 * t35, t14 (t14 * t29 + t7 * t25) * r_i_i_C(1) + (-t14 * t25 + t7 * t29) * r_i_i_C(2); 0 (t33 * t27 + t34 * t31) * t23, -t43, t47 * t12 + (-t41 * t26 - t30 * t43) * pkin(4) + t36 * (-t41 * t20 - t21 * t43) t45 (-t12 * t25 + t29 * t45) * r_i_i_C(1) + (-t12 * t29 - t25 * t45) * r_i_i_C(2);];
Ja_transl  = t5;
