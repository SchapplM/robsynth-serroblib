% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:50
% EndTime: 2019-02-26 19:56:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (110->34), mult. (286->63), div. (0->0), fcn. (362->10), ass. (0->28)
t32 = pkin(9) + r_i_i_C(3);
t14 = sin(pkin(6));
t18 = sin(qJ(4));
t31 = t14 * t18;
t19 = sin(qJ(2));
t30 = t14 * t19;
t21 = cos(qJ(4));
t29 = t14 * t21;
t22 = cos(qJ(2));
t28 = t14 * t22;
t16 = cos(pkin(6));
t27 = t16 * t19;
t26 = t16 * t22;
t17 = sin(qJ(5));
t20 = cos(qJ(5));
t25 = t20 * r_i_i_C(1) - t17 * r_i_i_C(2) + pkin(4);
t24 = t17 * r_i_i_C(1) + t20 * r_i_i_C(2) + pkin(2) + pkin(8);
t23 = t25 * t18 - t32 * t21 + qJ(3);
t15 = cos(pkin(11));
t13 = sin(pkin(11));
t11 = t16 * t21 - t18 * t28;
t9 = -t13 * t27 + t15 * t22;
t8 = t13 * t26 + t15 * t19;
t7 = t13 * t22 + t15 * t27;
t6 = t13 * t19 - t15 * t26;
t4 = t15 * t29 - t6 * t18;
t2 = t13 * t29 + t8 * t18;
t1 = [0, t23 * t9 - t24 * t8, t8, t32 * t2 + t25 * (-t13 * t31 + t8 * t21) (-t2 * t17 + t9 * t20) * r_i_i_C(1) + (-t9 * t17 - t2 * t20) * r_i_i_C(2), 0; 0, t23 * t7 - t24 * t6, t6, -t32 * t4 + t25 * (t15 * t31 + t6 * t21) (t4 * t17 + t7 * t20) * r_i_i_C(1) + (-t7 * t17 + t4 * t20) * r_i_i_C(2), 0; 1 (t23 * t19 + t24 * t22) * t14, -t28, t32 * t11 + t25 * (-t16 * t18 - t21 * t28) (-t11 * t17 + t20 * t30) * r_i_i_C(1) + (-t11 * t20 - t17 * t30) * r_i_i_C(2), 0;];
Ja_transl  = t1;
