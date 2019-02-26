% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRPRR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:47
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (60->34), mult. (173->69), div. (0->0), fcn. (226->12), ass. (0->28)
t13 = sin(pkin(11));
t15 = sin(pkin(6));
t27 = t13 * t15;
t20 = cos(pkin(6));
t26 = t13 * t20;
t14 = sin(pkin(7));
t25 = t14 * t15;
t18 = cos(pkin(11));
t24 = t18 * t15;
t23 = t18 * t20;
t11 = sin(pkin(13));
t16 = cos(pkin(13));
t21 = sin(qJ(3));
t22 = cos(qJ(3));
t10 = -t11 * t22 - t21 * t16;
t9 = t21 * t11 - t16 * t22;
t19 = cos(pkin(7));
t17 = cos(pkin(12));
t12 = sin(pkin(12));
t8 = -t12 * t26 + t17 * t18;
t7 = -t12 * t18 - t17 * t26;
t6 = t12 * t23 + t13 * t17;
t5 = -t12 * t13 + t17 * t23;
t4 = t10 * t19;
t3 = t9 * t19;
t2 = t10 * t14;
t1 = t9 * t14;
t28 = [0, t27 (-t1 * t27 + t10 * t8 - t3 * t7) * r_i_i_C(1) + (t2 * t27 + t4 * t7 + t8 * t9) * r_i_i_C(2) + (-t8 * t21 + (t13 * t25 + t19 * t7) * t22) * pkin(3), -t14 * t7 + t19 * t27, 0, 0; 0, -t24 (t1 * t24 + t10 * t6 - t3 * t5) * r_i_i_C(1) + (-t2 * t24 + t4 * t5 + t6 * t9) * r_i_i_C(2) + (-t6 * t21 + (-t14 * t24 + t19 * t5) * t22) * pkin(3), -t14 * t5 - t19 * t24, 0, 0; 1, t20 (pkin(3) * t14 * t22 - r_i_i_C(1) * t1 + r_i_i_C(2) * t2) * t20 + ((t10 * t12 - t17 * t3) * r_i_i_C(1) + (t12 * t9 + t17 * t4) * r_i_i_C(2) + (t17 * t19 * t22 - t12 * t21) * pkin(3)) * t15, -t17 * t25 + t19 * t20, 0, 0;];
Ja_transl  = t28;
