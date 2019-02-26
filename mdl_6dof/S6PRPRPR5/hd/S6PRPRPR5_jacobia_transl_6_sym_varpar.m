% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:45
% DurationCPUTime: 0.15s
% Computational Cost: add. (206->37), mult. (335->64), div. (0->0), fcn. (426->11), ass. (0->29)
t17 = pkin(11) + qJ(4);
t15 = sin(t17);
t16 = cos(t17);
t21 = sin(qJ(6));
t23 = cos(qJ(6));
t27 = t21 * r_i_i_C(1) + t23 * r_i_i_C(2) + qJ(5);
t31 = pkin(4) + pkin(9) + r_i_i_C(3);
t37 = t27 * t15 + t31 * t16 + cos(pkin(11)) * pkin(3) + pkin(2);
t18 = sin(pkin(10));
t19 = sin(pkin(6));
t36 = t18 * t19;
t22 = sin(qJ(2));
t35 = t19 * t22;
t24 = cos(qJ(2));
t34 = t19 * t24;
t33 = cos(pkin(6));
t32 = cos(pkin(10));
t30 = t18 * t33;
t29 = t19 * t32;
t28 = t33 * t32;
t26 = t23 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(5) + pkin(8) + qJ(3);
t10 = -t22 * t30 + t32 * t24;
t9 = t32 * t22 + t24 * t30;
t8 = t18 * t24 + t22 * t28;
t7 = t18 * t22 - t24 * t28;
t5 = t15 * t35 - t33 * t16;
t3 = t10 * t15 - t16 * t36;
t1 = t8 * t15 + t16 * t29;
t2 = [0, t26 * t10 - t37 * t9, t9, t27 * (t10 * t16 + t15 * t36) - t31 * t3, t3 (-t9 * t21 + t3 * t23) * r_i_i_C(1) + (-t3 * t21 - t9 * t23) * r_i_i_C(2); 0, t26 * t8 - t37 * t7, t7, t27 * (-t15 * t29 + t8 * t16) - t31 * t1, t1 (t1 * t23 - t7 * t21) * r_i_i_C(1) + (-t1 * t21 - t7 * t23) * r_i_i_C(2); 1 (t26 * t22 + t37 * t24) * t19, -t34, t27 * (t33 * t15 + t16 * t35) - t31 * t5, t5 (t21 * t34 + t5 * t23) * r_i_i_C(1) + (-t5 * t21 + t23 * t34) * r_i_i_C(2);];
Ja_transl  = t2;
