% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:57
% EndTime: 2019-02-26 21:56:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (48->21), mult. (109->31), div. (0->0), fcn. (128->6), ass. (0->23)
t12 = cos(qJ(2));
t24 = pkin(2) + pkin(3);
t9 = sin(qJ(2));
t14 = t9 * qJ(3) + t24 * t12;
t25 = pkin(1) + t14;
t8 = sin(qJ(4));
t23 = t12 * t8;
t13 = cos(qJ(1));
t22 = t13 * t9;
t20 = pkin(7) - pkin(8) - r_i_i_C(3);
t10 = sin(qJ(1));
t11 = cos(qJ(4));
t6 = t11 * t9 - t23;
t1 = t6 * t10;
t16 = t11 * t12 + t8 * t9;
t2 = t16 * t10;
t19 = t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
t3 = -t11 * t22 + t13 * t23;
t4 = t16 * t13;
t18 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
t17 = -t16 * r_i_i_C(1) - t6 * r_i_i_C(2);
t15 = qJ(3) * t12 - t24 * t9;
t5 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t25 * t10 + t20 * t13, t15 * t13 - t18, t22, t18, 0, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t20 * t10 + t25 * t13, t15 * t10 - t19, t10 * t9, t19, 0, 0; 0, t14 - t17, -t12, t17, 0, 0;];
Ja_transl  = t5;
