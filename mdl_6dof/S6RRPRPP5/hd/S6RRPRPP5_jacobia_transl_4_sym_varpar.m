% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_transl = S6RRPRPP5_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobia_transl_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:06
% EndTime: 2019-02-26 21:37:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (43->22), mult. (96->34), div. (0->0), fcn. (107->6), ass. (0->21)
t15 = pkin(2) + pkin(8) + r_i_i_C(3);
t9 = cos(qJ(2));
t13 = t15 * t9;
t6 = sin(qJ(2));
t21 = t6 * qJ(3) + pkin(1) + t13;
t7 = sin(qJ(1));
t20 = t7 * t6;
t8 = cos(qJ(4));
t19 = t7 * t8;
t18 = pkin(3) + pkin(7);
t10 = cos(qJ(1));
t17 = t10 * t6;
t16 = t10 * t8;
t5 = sin(qJ(4));
t12 = r_i_i_C(1) * t5 + r_i_i_C(2) * t8 + qJ(3);
t11 = t12 * t9 - t15 * t6;
t4 = -t5 * t20 + t16;
t3 = t10 * t5 + t6 * t19;
t2 = t5 * t17 + t19;
t1 = t6 * t16 - t7 * t5;
t14 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t18 * t10 - t21 * t7, t11 * t10, t17, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t21 * t10 + t18 * t7, t11 * t7, t20, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0, 0; 0, t12 * t6 + t13, -t9 (-r_i_i_C(1) * t8 + r_i_i_C(2) * t5) * t9, 0, 0;];
Ja_transl  = t14;
