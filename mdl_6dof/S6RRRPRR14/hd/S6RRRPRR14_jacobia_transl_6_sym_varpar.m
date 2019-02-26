% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR14_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR14_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:39
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.21s
% Computational Cost: add. (301->54), mult. (613->89), div. (0->0), fcn. (775->12), ass. (0->40)
t30 = sin(qJ(3));
t34 = cos(qJ(3));
t27 = qJ(5) + qJ(6);
t25 = sin(t27);
t26 = cos(t27);
t29 = sin(qJ(5));
t43 = pkin(5) * t29 + qJ(4);
t38 = r_i_i_C(1) * t25 + r_i_i_C(2) * t26 + t43;
t45 = pkin(3) + r_i_i_C(3) + pkin(11) + pkin(10);
t55 = t38 * t30 + t45 * t34 + pkin(2);
t31 = sin(qJ(2));
t32 = sin(qJ(1));
t35 = cos(qJ(2));
t46 = cos(pkin(6));
t50 = cos(qJ(1));
t41 = t46 * t50;
t17 = t31 * t32 - t35 * t41;
t18 = t31 * t41 + t32 * t35;
t28 = sin(pkin(6));
t44 = t28 * t50;
t9 = t18 * t30 + t34 * t44;
t54 = (-t17 * t25 + t26 * t9) * r_i_i_C(1) + (-t17 * t26 - t25 * t9) * r_i_i_C(2);
t42 = t32 * t46;
t20 = -t31 * t42 + t50 * t35;
t48 = t28 * t34;
t13 = t20 * t30 - t32 * t48;
t19 = t50 * t31 + t35 * t42;
t5 = t13 * t26 - t19 * t25;
t6 = t13 * t25 + t19 * t26;
t53 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t15 = t28 * t30 * t31 - t46 * t34;
t47 = t28 * t35;
t52 = (t15 * t26 + t25 * t47) * r_i_i_C(1) + (-t15 * t25 + t26 * t47) * r_i_i_C(2);
t33 = cos(qJ(5));
t51 = pkin(5) * t33 + pkin(4) + pkin(9);
t49 = t28 * t32;
t40 = r_i_i_C(1) * t26 - r_i_i_C(2) * t25 + t51;
t39 = t18 * t34 - t30 * t44;
t14 = t20 * t34 + t30 * t49;
t1 = [-t32 * pkin(1) - t18 * pkin(2) + pkin(8) * t44 - t40 * t17 - t38 * t9 - t45 * t39, -t19 * t55 + t40 * t20, -t45 * t13 + t38 * t14, t13 (t13 * t33 - t19 * t29) * pkin(5) + t53, t53; t50 * pkin(1) + t20 * pkin(2) + pkin(8) * t49 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t43 * t13 + t45 * t14 + t51 * t19, -t17 * t55 + t40 * t18, t38 * t39 - t45 * t9, t9 (-t17 * t29 + t33 * t9) * pkin(5) + t54, t54; 0 (t40 * t31 + t55 * t35) * t28, -t45 * t15 + t38 * (t46 * t30 + t31 * t48) t15 (t15 * t33 + t29 * t47) * pkin(5) + t52, t52;];
Ja_transl  = t1;
