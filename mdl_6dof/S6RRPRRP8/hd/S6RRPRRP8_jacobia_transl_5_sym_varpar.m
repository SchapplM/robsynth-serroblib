% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:52
% EndTime: 2019-02-26 21:49:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (152->32), mult. (133->47), div. (0->0), fcn. (148->10), ass. (0->28)
t20 = cos(qJ(2));
t18 = sin(qJ(2));
t29 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
t25 = t29 * t18;
t17 = pkin(10) + qJ(4);
t14 = cos(t17);
t9 = pkin(4) * t14 + cos(pkin(10)) * pkin(3) + pkin(2);
t34 = t20 * t9 + pkin(1) + t25;
t15 = qJ(5) + t17;
t11 = sin(t15);
t12 = cos(t15);
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t28 = t19 * t20;
t5 = t11 * t28 + t12 * t21;
t6 = t11 * t21 - t12 * t28;
t33 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t27 = t20 * t21;
t7 = -t11 * t27 + t19 * t12;
t8 = t19 * t11 + t12 * t27;
t32 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t13 = sin(t17);
t31 = pkin(4) * t13;
t30 = pkin(7) + t31 + sin(pkin(10)) * pkin(3);
t24 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
t23 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
t22 = -t23 * t18 + t29 * t20;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t34 * t19 + t30 * t21, t22 * t21, t21 * t18 (-t13 * t27 + t14 * t19) * pkin(4) + t32, t32, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t30 * t19 + t34 * t21, t22 * t19, t19 * t18 (-t13 * t28 - t14 * t21) * pkin(4) + t33, t33, 0; 0, t23 * t20 + t25, -t20 (t24 - t31) * t18, t24 * t18, 0;];
Ja_transl  = t1;
