% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR6
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

function Ja_transl = S6RPRRRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:36
% EndTime: 2019-02-26 21:17:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (247->39), mult. (170->51), div. (0->0), fcn. (188->11), ass. (0->35)
t22 = pkin(11) + qJ(3);
t17 = cos(t22);
t16 = sin(t22);
t40 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t32 = t40 * t16;
t24 = qJ(4) + qJ(5);
t19 = cos(t24);
t11 = pkin(5) * t19 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t44 = t32 + t17 * t9 + cos(pkin(11)) * pkin(2) + pkin(1);
t20 = qJ(6) + t24;
t15 = cos(t20);
t27 = cos(qJ(1));
t34 = t27 * t15;
t14 = sin(t20);
t26 = sin(qJ(1));
t37 = t26 * t14;
t5 = t17 * t37 + t34;
t35 = t27 * t14;
t36 = t26 * t15;
t6 = -t17 * t36 + t35;
t43 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t17 * t35 + t36;
t8 = t17 * t34 + t37;
t42 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t18 = sin(t24);
t41 = pkin(5) * t18;
t39 = t17 * t26;
t38 = t17 * t27;
t10 = t41 + sin(qJ(4)) * pkin(4);
t33 = t10 + pkin(7) + qJ(2);
t30 = -r_i_i_C(1) * t14 - r_i_i_C(2) * t15;
t29 = r_i_i_C(1) * t15 - r_i_i_C(2) * t14 + t9;
t28 = -t29 * t16 + t40 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t44 * t26 + t33 * t27, t26, t28 * t27, -t10 * t38 + t26 * t11 + t42 (-t18 * t38 + t19 * t26) * pkin(5) + t42, t42; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t33 * t26 + t44 * t27, -t27, t28 * t26, -t10 * t39 - t27 * t11 + t43 (-t18 * t39 - t19 * t27) * pkin(5) + t43, t43; 0, 0, t29 * t17 + t32 (-t10 + t30) * t16 (t30 - t41) * t16, t30 * t16;];
Ja_transl  = t1;
