% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:59
% EndTime: 2019-02-26 20:00:00
% DurationCPUTime: 0.22s
% Computational Cost: add. (246->53), mult. (658->98), div. (0->0), fcn. (853->12), ass. (0->43)
t31 = sin(qJ(3));
t47 = -r_i_i_C(3) - pkin(9) + qJ(4);
t53 = cos(qJ(3));
t54 = pkin(3) * t53 + t31 * t47 + pkin(2);
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t52 = t28 * t32;
t34 = cos(qJ(2));
t51 = t28 * t34;
t50 = cos(pkin(6));
t49 = cos(pkin(10));
t48 = sin(pkin(10));
t27 = sin(pkin(11));
t46 = t27 * t53;
t29 = cos(pkin(11));
t45 = t29 * t53;
t44 = t34 * t53;
t43 = t28 * t49;
t42 = t28 * t48;
t40 = t50 * t49;
t39 = t50 * t48;
t30 = sin(qJ(6));
t33 = cos(qJ(6));
t38 = r_i_i_C(1) * t30 + r_i_i_C(2) * t33 + qJ(5);
t37 = r_i_i_C(1) * t33 - r_i_i_C(2) * t30 + pkin(4) + pkin(5);
t35 = -t27 * t38 - t29 * t37 - pkin(3);
t22 = t31 * t50 + t52 * t53;
t21 = t31 * t52 - t50 * t53;
t20 = -t32 * t39 + t34 * t49;
t19 = t32 * t49 + t34 * t39;
t18 = t32 * t40 + t34 * t48;
t17 = t32 * t48 - t34 * t40;
t14 = t20 * t53 + t31 * t42;
t13 = t20 * t31 - t42 * t53;
t12 = t18 * t53 - t31 * t43;
t11 = t18 * t31 + t43 * t53;
t10 = t22 * t29 - t27 * t51;
t9 = t22 * t27 + t29 * t51;
t4 = t14 * t29 + t19 * t27;
t3 = t14 * t27 - t19 * t29;
t2 = t12 * t29 + t17 * t27;
t1 = t12 * t27 - t17 * t29;
t5 = [0, t20 * pkin(8) + t38 * (-t19 * t46 - t20 * t29) + t37 * (-t19 * t45 + t20 * t27) - t54 * t19, t13 * t35 + t14 * t47, t13, t3 (t3 * t33 - t30 * t4) * r_i_i_C(1) + (-t3 * t30 - t33 * t4) * r_i_i_C(2); 0, t18 * pkin(8) + t38 * (-t17 * t46 - t18 * t29) + t37 * (-t17 * t45 + t18 * t27) - t54 * t17, t11 * t35 + t12 * t47, t11, t1 (t1 * t33 - t2 * t30) * r_i_i_C(1) + (-t1 * t30 - t2 * t33) * r_i_i_C(2); 1 (t38 * (t27 * t44 - t29 * t32) + t37 * (t27 * t32 + t29 * t44) + t32 * pkin(8) + t54 * t34) * t28, t21 * t35 + t22 * t47, t21, t9 (-t10 * t30 + t33 * t9) * r_i_i_C(1) + (-t10 * t33 - t30 * t9) * r_i_i_C(2);];
Ja_transl  = t5;
