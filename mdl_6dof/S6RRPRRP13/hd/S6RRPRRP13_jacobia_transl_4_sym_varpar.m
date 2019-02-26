% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP13_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP13_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:58
% EndTime: 2019-02-26 21:52:58
% DurationCPUTime: 0.09s
% Computational Cost: add. (75->31), mult. (179->50), div. (0->0), fcn. (222->8), ass. (0->26)
t27 = pkin(3) + pkin(8);
t10 = sin(pkin(6));
t14 = sin(qJ(1));
t26 = t10 * t14;
t16 = cos(qJ(2));
t25 = t10 * t16;
t17 = cos(qJ(1));
t24 = t10 * t17;
t13 = sin(qJ(2));
t23 = t13 * t14;
t22 = t13 * t17;
t21 = t14 * t16;
t20 = t16 * t17;
t19 = -r_i_i_C(3) - pkin(9) - pkin(2);
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t18 = t12 * r_i_i_C(1) + t15 * r_i_i_C(2) + qJ(3);
t11 = cos(pkin(6));
t8 = t15 * t24;
t6 = -t11 * t23 + t20;
t5 = t11 * t21 + t22;
t4 = t11 * t22 + t21;
t3 = -t11 * t20 + t23;
t2 = t12 * t5 + t15 * t26;
t1 = -t12 * t26 + t15 * t5;
t7 = [-t14 * pkin(1) + t8 * r_i_i_C(1) - t18 * t3 + (-t12 * r_i_i_C(2) + t27) * t24 + t19 * t4, t18 * t6 + t19 * t5, t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2), 0, 0; pkin(1) * t17 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(3) - t19 * t6 + t27 * t26, t18 * t4 + t19 * t3, t3 (t12 * t24 + t3 * t15) * r_i_i_C(1) + (-t3 * t12 + t8) * r_i_i_C(2), 0, 0; 0 (t18 * t13 - t19 * t16) * t10, -t25 (-t11 * t12 - t15 * t25) * r_i_i_C(1) + (-t11 * t15 + t12 * t25) * r_i_i_C(2), 0, 0;];
Ja_transl  = t7;
