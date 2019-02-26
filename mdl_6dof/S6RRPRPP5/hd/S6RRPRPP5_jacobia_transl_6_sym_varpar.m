% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:06
% EndTime: 2019-02-26 21:37:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (89->28), mult. (192->37), div. (0->0), fcn. (223->6), ass. (0->23)
t10 = cos(qJ(2));
t15 = pkin(2) + pkin(8) - r_i_i_C(3) - qJ(6);
t14 = t15 * t10;
t7 = sin(qJ(2));
t24 = t7 * qJ(3) + pkin(1) + t14;
t8 = sin(qJ(1));
t23 = t8 * t7;
t9 = cos(qJ(4));
t22 = t8 * t9;
t21 = pkin(3) + pkin(7);
t11 = cos(qJ(1));
t20 = t11 * t7;
t19 = t11 * t9;
t18 = r_i_i_C(2) + qJ(5);
t17 = pkin(4) + pkin(5) + r_i_i_C(1);
t6 = sin(qJ(4));
t13 = t17 * t6 - t18 * t9 + qJ(3);
t12 = t13 * t10 - t15 * t7;
t4 = -t6 * t23 + t19;
t3 = t11 * t6 + t7 * t22;
t2 = t6 * t20 + t22;
t1 = -t7 * t19 + t8 * t6;
t5 = [t21 * t11 + t17 * t4 + t18 * t3 - t24 * t8, t12 * t11, t20, -t17 * t1 + t18 * t2, t1, -t11 * t10; t18 * t1 + t24 * t11 + t17 * t2 + t21 * t8, t12 * t8, t23, t17 * t3 - t18 * t4, -t3, -t8 * t10; 0, t13 * t7 + t14, -t10 (-t17 * t9 - t18 * t6) * t10, t10 * t9, -t7;];
Ja_transl  = t5;
