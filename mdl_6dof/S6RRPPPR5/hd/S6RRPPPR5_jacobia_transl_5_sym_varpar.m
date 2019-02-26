% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:12
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (59->20), mult. (136->31), div. (0->0), fcn. (159->6), ass. (0->20)
t10 = cos(qJ(2));
t15 = pkin(4) - r_i_i_C(2) + qJ(3);
t8 = sin(qJ(2));
t13 = t15 * t8;
t21 = t10 * pkin(2) + pkin(1) + t13;
t16 = pkin(3) + r_i_i_C(3) + qJ(5);
t17 = r_i_i_C(1) + qJ(4);
t6 = sin(pkin(9));
t7 = cos(pkin(9));
t20 = t16 * t7 + t17 * t6 + pkin(2);
t9 = sin(qJ(1));
t19 = t10 * t9;
t11 = cos(qJ(1));
t18 = t10 * t11;
t12 = t15 * t10 - t20 * t8;
t4 = t7 * t18 + t9 * t6;
t3 = t6 * t18 - t9 * t7;
t2 = -t11 * t6 + t7 * t19;
t1 = t11 * t7 + t6 * t19;
t5 = [pkin(7) * t11 - t17 * t1 - t16 * t2 - t21 * t9, t12 * t11, t11 * t8, t3, t4, 0; t9 * pkin(7) + t21 * t11 + t16 * t4 + t17 * t3, t12 * t9, t9 * t8, t1, t2, 0; 0, t20 * t10 + t13, -t10, t8 * t6, t8 * t7, 0;];
Ja_transl  = t5;
