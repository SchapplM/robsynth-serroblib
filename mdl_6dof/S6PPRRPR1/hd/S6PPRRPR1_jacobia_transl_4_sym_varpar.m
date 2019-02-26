% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:18
% EndTime: 2019-02-26 19:40:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (102->33), mult. (288->67), div. (0->0), fcn. (377->12), ass. (0->34)
t36 = pkin(9) + r_i_i_C(3);
t15 = sin(pkin(11));
t17 = sin(pkin(6));
t35 = t15 * t17;
t21 = cos(pkin(6));
t34 = t15 * t21;
t16 = sin(pkin(7));
t33 = t16 * t17;
t32 = t16 * t21;
t18 = cos(pkin(12));
t20 = cos(pkin(7));
t31 = t18 * t20;
t19 = cos(pkin(11));
t30 = t19 * t17;
t29 = t19 * t21;
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t28 = t24 * r_i_i_C(1) - t22 * r_i_i_C(2) + pkin(3);
t14 = sin(pkin(12));
t10 = -t15 * t14 + t18 * t29;
t27 = t10 * t20 - t16 * t30;
t12 = -t19 * t14 - t18 * t34;
t26 = t12 * t20 + t15 * t33;
t25 = cos(qJ(3));
t23 = sin(qJ(3));
t13 = -t14 * t34 + t19 * t18;
t11 = t14 * t29 + t15 * t18;
t9 = -t18 * t33 + t21 * t20;
t8 = -t12 * t16 + t20 * t35;
t7 = -t10 * t16 - t20 * t30;
t6 = t23 * t32 + (t14 * t25 + t23 * t31) * t17;
t4 = t13 * t25 + t26 * t23;
t2 = t11 * t25 + t27 * t23;
t1 = [0, t35, t36 * t4 + t28 * (-t13 * t23 + t26 * t25) (-t4 * t22 + t8 * t24) * r_i_i_C(1) + (-t8 * t22 - t4 * t24) * r_i_i_C(2), 0, 0; 0, -t30, t36 * t2 + t28 * (-t11 * t23 + t27 * t25) (-t2 * t22 + t7 * t24) * r_i_i_C(1) + (-t2 * t24 - t7 * t22) * r_i_i_C(2), 0, 0; 1, t21, t36 * t6 + t28 * (t25 * t32 + (-t14 * t23 + t25 * t31) * t17) (-t6 * t22 + t9 * t24) * r_i_i_C(1) + (-t9 * t22 - t6 * t24) * r_i_i_C(2), 0, 0;];
Ja_transl  = t1;
