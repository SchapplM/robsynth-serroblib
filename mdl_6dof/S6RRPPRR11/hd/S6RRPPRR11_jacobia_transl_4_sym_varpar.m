% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR11_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->20), mult. (150->32), div. (0->0), fcn. (188->8), ass. (0->18)
t19 = cos(pkin(6));
t18 = r_i_i_C(3) + qJ(4) + pkin(2);
t11 = sin(qJ(1));
t17 = t11 * t19;
t13 = cos(qJ(1));
t16 = t13 * t19;
t7 = sin(pkin(11));
t9 = cos(pkin(11));
t15 = r_i_i_C(1) * t7 + r_i_i_C(2) * t9 + qJ(3);
t8 = sin(pkin(6));
t14 = (r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + pkin(3) + pkin(8)) * t8;
t12 = cos(qJ(2));
t10 = sin(qJ(2));
t4 = -t10 * t17 + t12 * t13;
t3 = t13 * t10 + t12 * t17;
t2 = t10 * t16 + t11 * t12;
t1 = t10 * t11 - t12 * t16;
t5 = [-t11 * pkin(1) - t15 * t1 + t13 * t14 - t18 * t2, t15 * t4 - t18 * t3, t3, t4, 0, 0; t13 * pkin(1) + t11 * t14 + t15 * t3 + t18 * t4, -t18 * t1 + t15 * t2, t1, t2, 0, 0; 0 (t15 * t10 + t18 * t12) * t8, -t8 * t12, t8 * t10, 0, 0;];
Ja_transl  = t5;
