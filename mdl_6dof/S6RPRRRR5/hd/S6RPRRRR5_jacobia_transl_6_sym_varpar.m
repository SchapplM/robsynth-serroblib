% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:10
% EndTime: 2019-02-26 21:17:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (244->42), mult. (167->53), div. (0->0), fcn. (181->11), ass. (0->40)
t26 = pkin(11) + qJ(3);
t22 = qJ(4) + t26;
t18 = sin(t22);
t19 = cos(t22);
t30 = cos(qJ(5));
t20 = t30 * pkin(5) + pkin(4);
t32 = -pkin(10) - pkin(9);
t55 = t19 * t20 + (r_i_i_C(3) - t32) * t18;
t17 = pkin(3) * cos(t26);
t54 = t17 + cos(pkin(11)) * pkin(2) + pkin(1) + t55;
t27 = qJ(5) + qJ(6);
t24 = cos(t27);
t31 = cos(qJ(1));
t41 = t31 * t24;
t23 = sin(t27);
t29 = sin(qJ(1));
t44 = t29 * t23;
t5 = t19 * t44 + t41;
t42 = t31 * t23;
t43 = t29 * t24;
t6 = -t19 * t43 + t42;
t53 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t19 * t42 + t43;
t8 = t19 * t41 + t44;
t52 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t28 = sin(qJ(5));
t51 = pkin(5) * t28;
t50 = r_i_i_C(1) * t24;
t49 = r_i_i_C(2) * t23;
t39 = t18 * t49;
t46 = t19 * t29;
t47 = r_i_i_C(3) * t46 + t29 * t39;
t45 = t19 * t31;
t40 = r_i_i_C(3) * t45 + t31 * t39;
t38 = pkin(8) + pkin(7) + qJ(2) + t51;
t36 = -r_i_i_C(1) * t23 - r_i_i_C(2) * t24;
t35 = -t19 * t32 + (-t20 - t50) * t18;
t34 = (-t49 + t50) * t19 + t55;
t33 = -pkin(3) * sin(t26) + t35;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t54 * t29 + t38 * t31, t29, t31 * t33 + t40, t31 * t35 + t40 (-t28 * t45 + t29 * t30) * pkin(5) + t52, t52; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t38 * t29 + t54 * t31, -t31, t29 * t33 + t47, t29 * t35 + t47 (-t28 * t46 - t30 * t31) * pkin(5) + t53, t53; 0, 0, t17 + t34, t34 (t36 - t51) * t18, t36 * t18;];
Ja_transl  = t1;
