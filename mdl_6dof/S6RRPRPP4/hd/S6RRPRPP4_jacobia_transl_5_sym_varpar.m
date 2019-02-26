% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:40
% EndTime: 2019-02-26 21:36:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (84->29), mult. (126->44), div. (0->0), fcn. (140->8), ass. (0->24)
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t20 = pkin(2) + r_i_i_C(3) + qJ(5) + pkin(8);
t18 = t20 * t14;
t10 = sin(qJ(4));
t19 = pkin(4) * t10 + qJ(3);
t26 = -t19 * t11 - pkin(1) - t18;
t13 = cos(qJ(4));
t23 = pkin(4) * t13;
t24 = pkin(7) + pkin(3) + t23;
t12 = sin(qJ(1));
t22 = t12 * t11;
t15 = cos(qJ(1));
t21 = t15 * t11;
t8 = qJ(4) + pkin(9);
t6 = sin(t8);
t7 = cos(t8);
t17 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + t19;
t16 = -t20 * t11 + t17 * t14;
t4 = t15 * t7 - t6 * t22;
t3 = t15 * t6 + t7 * t22;
t2 = t12 * t7 + t6 * t21;
t1 = -t12 * t6 + t7 * t21;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t26 * t12 + t24 * t15, t16 * t15, t21, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (-t10 * t12 + t13 * t21) * pkin(4), t15 * t14, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t24 * t12 - t26 * t15, t16 * t12, t22, t3 * r_i_i_C(1) + t4 * r_i_i_C(2) + (t10 * t15 + t13 * t22) * pkin(4), t12 * t14, 0; 0, t17 * t11 + t18, -t14 (-r_i_i_C(1) * t7 + r_i_i_C(2) * t6 - t23) * t14, t11, 0;];
Ja_transl  = t5;
