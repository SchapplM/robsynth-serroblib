% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:39
% EndTime: 2019-02-26 21:18:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (182->42), mult. (169->54), div. (0->0), fcn. (183->10), ass. (0->38)
t22 = qJ(3) + qJ(4);
t20 = cos(t22);
t29 = -pkin(10) - pkin(9);
t55 = (r_i_i_C(3) - t29) * t20;
t18 = sin(t22);
t21 = qJ(5) + qJ(6);
t17 = sin(t21);
t48 = r_i_i_C(2) * t17;
t53 = t18 * t29 + t20 * t48;
t19 = cos(t21);
t28 = cos(qJ(1));
t41 = t28 * t19;
t25 = sin(qJ(1));
t44 = t25 * t17;
t5 = -t18 * t44 + t41;
t42 = t28 * t17;
t43 = t25 * t19;
t6 = t18 * t43 + t42;
t52 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t7 = t18 * t42 + t43;
t8 = t18 * t41 - t44;
t51 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t50 = cos(qJ(3)) * pkin(3);
t23 = sin(qJ(5));
t49 = pkin(5) * t23;
t47 = sin(qJ(3)) * pkin(3);
t46 = t18 * t25;
t40 = t53 * t28;
t26 = cos(qJ(5));
t16 = t26 * pkin(5) + pkin(4);
t38 = r_i_i_C(3) * t46 + (r_i_i_C(1) * t43 + t16 * t25) * t20;
t37 = -r_i_i_C(1) * t19 - t16;
t36 = pkin(1) + pkin(8) + pkin(7) + t49;
t35 = -r_i_i_C(1) * t17 - r_i_i_C(2) * t19;
t33 = -r_i_i_C(3) * t18 + t37 * t20;
t32 = t55 + (t37 + t48) * t18;
t31 = t18 * t16 + qJ(2) + t47 - t55;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t36 * t25 + t31 * t28, t25 (-t53 + t50) * t25 + t38, -t25 * t53 + t38 (-t23 * t46 + t26 * t28) * pkin(5) + t52, t52; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t31 * t25 + t36 * t28, -t28 (t33 - t50) * t28 + t40, t33 * t28 + t40 (t18 * t23 * t28 + t25 * t26) * pkin(5) + t51, t51; 0, 0, t32 - t47, t32 (t35 - t49) * t20, t35 * t20;];
Ja_transl  = t1;
