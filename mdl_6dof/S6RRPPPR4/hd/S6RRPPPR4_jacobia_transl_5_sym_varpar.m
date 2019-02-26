% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:35
% EndTime: 2019-02-26 21:23:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (52->22), mult. (117->30), div. (0->0), fcn. (136->6), ass. (0->19)
t15 = pkin(2) + r_i_i_C(2) + qJ(4);
t9 = cos(qJ(2));
t13 = t15 * t9;
t7 = sin(qJ(2));
t21 = t7 * qJ(3) + pkin(1) + t13;
t8 = sin(qJ(1));
t20 = t8 * t7;
t19 = pkin(3) + pkin(7);
t18 = pkin(4) + r_i_i_C(1);
t10 = cos(qJ(1));
t17 = t10 * t7;
t16 = r_i_i_C(3) + qJ(5);
t5 = sin(pkin(9));
t6 = cos(pkin(9));
t12 = -t16 * t6 + t18 * t5 + qJ(3);
t11 = t12 * t9 - t15 * t7;
t3 = t10 * t5 + t6 * t20;
t1 = -t6 * t17 + t8 * t5;
t2 = [t18 * (t10 * t6 - t5 * t20) + t16 * t3 + t19 * t10 - t21 * t8, t11 * t10, t17, t10 * t9, t1, 0; t19 * t8 + t18 * (t5 * t17 + t8 * t6) + t16 * t1 + t21 * t10, t11 * t8, t20, t8 * t9, -t3, 0; 0, t12 * t7 + t13, -t9, t7, t9 * t6, 0;];
Ja_transl  = t2;
