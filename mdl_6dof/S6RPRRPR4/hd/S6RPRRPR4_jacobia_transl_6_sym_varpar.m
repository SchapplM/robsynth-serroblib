% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:37
% EndTime: 2019-02-26 21:02:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (211->39), mult. (137->45), div. (0->0), fcn. (150->11), ass. (0->32)
t23 = pkin(10) + qJ(3);
t20 = qJ(4) + t23;
t14 = sin(t20);
t15 = cos(t20);
t22 = pkin(11) + qJ(6);
t17 = sin(t22);
t42 = r_i_i_C(2) * t17;
t46 = r_i_i_C(3) * t15 + t14 * t42;
t16 = cos(pkin(11)) * pkin(5) + pkin(4);
t25 = -pkin(9) - qJ(5);
t45 = t15 * t16 + (r_i_i_C(3) - t25) * t14;
t13 = pkin(3) * cos(t23);
t44 = t13 + cos(pkin(10)) * pkin(2) + pkin(1) + t45;
t19 = cos(t22);
t43 = r_i_i_C(1) * t19;
t26 = sin(qJ(1));
t39 = t46 * t26;
t27 = cos(qJ(1));
t38 = t46 * t27;
t37 = t17 * t27;
t36 = t19 * t27;
t35 = t26 * t17;
t34 = t26 * t19;
t32 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(2);
t30 = -t15 * t25 + (-t16 - t43) * t14;
t29 = (-t42 + t43) * t15 + t45;
t28 = -pkin(3) * sin(t23) + t30;
t4 = t15 * t36 + t35;
t3 = -t15 * t37 + t34;
t2 = -t15 * t34 + t37;
t1 = t15 * t35 + t36;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t44 * t26 + t32 * t27, t26, t28 * t27 + t38, t30 * t27 + t38, t27 * t14, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t32 * t26 + t44 * t27, -t27, t28 * t26 + t39, t30 * t26 + t39, t26 * t14, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, t13 + t29, t29, -t15 (-r_i_i_C(1) * t17 - r_i_i_C(2) * t19) * t14;];
Ja_transl  = t5;
