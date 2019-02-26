% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:14
% EndTime: 2019-02-26 19:48:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (226->39), mult. (313->66), div. (0->0), fcn. (397->13), ass. (0->33)
t21 = pkin(11) + qJ(4);
t17 = sin(t21);
t19 = cos(t21);
t20 = pkin(12) + qJ(6);
t16 = sin(t20);
t18 = cos(t20);
t31 = t18 * r_i_i_C(1) - t16 * r_i_i_C(2) + cos(pkin(12)) * pkin(5) + pkin(4);
t40 = r_i_i_C(3) + pkin(9) + qJ(5);
t41 = t40 * t17 + t31 * t19 + cos(pkin(11)) * pkin(3) + pkin(2);
t23 = sin(pkin(10));
t24 = sin(pkin(6));
t39 = t23 * t24;
t27 = sin(qJ(2));
t38 = t24 * t27;
t28 = cos(qJ(2));
t37 = t24 * t28;
t36 = cos(pkin(6));
t35 = cos(pkin(10));
t34 = t23 * t36;
t33 = t24 * t35;
t32 = t36 * t35;
t30 = sin(pkin(12)) * pkin(5) + t16 * r_i_i_C(1) + t18 * r_i_i_C(2) + pkin(8) + qJ(3);
t10 = -t27 * t34 + t35 * t28;
t9 = t35 * t27 + t28 * t34;
t8 = t23 * t28 + t27 * t32;
t7 = t23 * t27 - t28 * t32;
t6 = t36 * t17 + t19 * t38;
t5 = t17 * t38 - t36 * t19;
t4 = t10 * t19 + t17 * t39;
t3 = t10 * t17 - t19 * t39;
t2 = -t17 * t33 + t8 * t19;
t1 = t8 * t17 + t19 * t33;
t11 = [0, t30 * t10 - t41 * t9, t9, -t31 * t3 + t40 * t4, t3 (-t4 * t16 + t9 * t18) * r_i_i_C(1) + (-t9 * t16 - t4 * t18) * r_i_i_C(2); 0, t30 * t8 - t41 * t7, t7, -t31 * t1 + t40 * t2, t1 (-t2 * t16 + t7 * t18) * r_i_i_C(1) + (-t7 * t16 - t2 * t18) * r_i_i_C(2); 1 (t30 * t27 + t41 * t28) * t24, -t37, -t31 * t5 + t40 * t6, t5 (-t6 * t16 - t18 * t37) * r_i_i_C(1) + (t16 * t37 - t6 * t18) * r_i_i_C(2);];
Ja_transl  = t11;
