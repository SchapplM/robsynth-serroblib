% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6RRPRPP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobia_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:06
% EndTime: 2019-02-26 21:37:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (69->25), mult. (155->35), div. (0->0), fcn. (179->6), ass. (0->23)
t10 = cos(qJ(2));
t16 = pkin(2) + pkin(8) + r_i_i_C(2);
t14 = t16 * t10;
t7 = sin(qJ(2));
t24 = qJ(3) * t7 + pkin(1) + t14;
t8 = sin(qJ(1));
t23 = t8 * t7;
t9 = cos(qJ(4));
t22 = t8 * t9;
t21 = pkin(3) + pkin(7);
t20 = pkin(4) + r_i_i_C(1);
t11 = cos(qJ(1));
t19 = t11 * t7;
t18 = t11 * t9;
t17 = r_i_i_C(3) + qJ(5);
t6 = sin(qJ(4));
t13 = -t17 * t9 + t20 * t6 + qJ(3);
t12 = t13 * t10 - t16 * t7;
t4 = -t6 * t23 + t18;
t3 = t11 * t6 + t7 * t22;
t2 = t6 * t19 + t22;
t1 = -t7 * t18 + t6 * t8;
t5 = [t21 * t11 + t17 * t3 + t20 * t4 - t24 * t8, t12 * t11, t19, -t20 * t1 + t17 * t2, t1, 0; t17 * t1 + t24 * t11 + t20 * t2 + t21 * t8, t12 * t8, t23, -t17 * t4 + t20 * t3, -t3, 0; 0, t13 * t7 + t14, -t10 (-t17 * t6 - t20 * t9) * t10, t10 * t9, 0;];
Ja_transl  = t5;
