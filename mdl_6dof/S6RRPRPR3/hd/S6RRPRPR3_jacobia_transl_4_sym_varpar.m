% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:00
% EndTime: 2019-02-26 21:39:00
% DurationCPUTime: 0.11s
% Computational Cost: add. (74->25), mult. (90->34), div. (0->0), fcn. (100->8), ass. (0->22)
t24 = pkin(8) + r_i_i_C(3);
t9 = qJ(2) + pkin(10);
t6 = sin(t9);
t26 = t24 * t6 + cos(qJ(2)) * pkin(2);
t7 = cos(t9);
t25 = t7 * pkin(3) + pkin(1) + t26;
t11 = sin(qJ(4));
t15 = cos(qJ(1));
t23 = t11 * t15;
t13 = sin(qJ(1));
t22 = t13 * t11;
t14 = cos(qJ(4));
t21 = t13 * t14;
t20 = t14 * t15;
t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t11 + pkin(3);
t16 = -sin(qJ(2)) * pkin(2) - t17 * t6 + t24 * t7;
t10 = -qJ(3) - pkin(7);
t4 = t7 * t20 + t22;
t3 = -t7 * t23 + t21;
t2 = -t7 * t21 + t23;
t1 = t7 * t22 + t20;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t15 * t10 - t25 * t13, t16 * t15, t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t13 * t10 + t25 * t15, t16 * t13, -t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t17 * t7 + t26, 0 (-r_i_i_C(1) * t11 - r_i_i_C(2) * t14) * t6, 0, 0;];
Ja_transl  = t5;
