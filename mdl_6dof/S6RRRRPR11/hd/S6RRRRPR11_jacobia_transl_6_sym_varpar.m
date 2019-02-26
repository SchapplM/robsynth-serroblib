% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:21
% DurationCPUTime: 0.21s
% Computational Cost: add. (362->56), mult. (542->88), div. (0->0), fcn. (682->14), ass. (0->42)
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t34 = qJ(4) + pkin(12);
t23 = pkin(5) * cos(t34) + cos(qJ(4)) * pkin(4);
t21 = pkin(3) + t23;
t31 = qJ(6) + t34;
t29 = sin(t31);
t30 = cos(t31);
t43 = r_i_i_C(1) * t30 - r_i_i_C(2) * t29 + t21;
t52 = r_i_i_C(3) + pkin(11) + qJ(5) + pkin(10);
t57 = t52 * t36 + t43 * t39 + pkin(2);
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t40 = cos(qJ(2));
t47 = cos(pkin(6));
t51 = cos(qJ(1));
t44 = t47 * t51;
t18 = t37 * t44 + t38 * t40;
t35 = sin(pkin(6));
t46 = t35 * t51;
t10 = t18 * t39 - t36 * t46;
t17 = t38 * t37 - t40 * t44;
t56 = (-t10 * t29 + t17 * t30) * r_i_i_C(1) + (-t10 * t30 - t17 * t29) * r_i_i_C(2);
t45 = t38 * t47;
t20 = -t37 * t45 + t51 * t40;
t50 = t35 * t38;
t14 = t20 * t39 + t36 * t50;
t19 = t51 * t37 + t40 * t45;
t5 = -t14 * t29 + t19 * t30;
t6 = t14 * t30 + t19 * t29;
t55 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t49 = t35 * t39;
t16 = t47 * t36 + t37 * t49;
t48 = t35 * t40;
t54 = (-t16 * t29 - t30 * t48) * r_i_i_C(1) + (-t16 * t30 + t29 * t48) * r_i_i_C(2);
t22 = pkin(5) * sin(t34) + sin(qJ(4)) * pkin(4);
t53 = pkin(9) + t22;
t42 = t29 * r_i_i_C(1) + t30 * r_i_i_C(2) + t53;
t9 = t18 * t36 + t39 * t46;
t15 = t35 * t37 * t36 - t47 * t39;
t13 = t20 * t36 - t38 * t49;
t1 = [-t38 * pkin(1) - t18 * pkin(2) + pkin(8) * t46 - t43 * t10 - t42 * t17 - t52 * t9, -t19 * t57 + t42 * t20, -t43 * t13 + t52 * t14, -t14 * t22 + t19 * t23 + t55, t13, t55; t51 * pkin(1) + t20 * pkin(2) + pkin(8) * t50 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t52 * t13 + t14 * t21 + t53 * t19, -t17 * t57 + t42 * t18, t52 * t10 - t43 * t9, -t10 * t22 + t17 * t23 + t56, t9, t56; 0 (t42 * t37 + t57 * t40) * t35, -t43 * t15 + t52 * t16, -t16 * t22 - t23 * t48 + t54, t15, t54;];
Ja_transl  = t1;
