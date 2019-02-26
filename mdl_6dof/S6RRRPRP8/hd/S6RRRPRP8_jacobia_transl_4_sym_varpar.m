% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP8_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP8_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobia_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:08
% EndTime: 2019-02-26 22:13:08
% DurationCPUTime: 0.12s
% Computational Cost: add. (61->22), mult. (142->34), div. (0->0), fcn. (163->6), ass. (0->21)
t11 = cos(qJ(2));
t19 = pkin(8) + r_i_i_C(2);
t8 = sin(qJ(2));
t15 = t19 * t8;
t23 = t11 * pkin(2) + pkin(1) + t15;
t10 = cos(qJ(3));
t16 = r_i_i_C(3) + qJ(4);
t20 = pkin(3) + r_i_i_C(1);
t7 = sin(qJ(3));
t22 = t20 * t10 + t16 * t7 + pkin(2);
t9 = sin(qJ(1));
t21 = t9 * t7;
t18 = t9 * t10;
t12 = cos(qJ(1));
t17 = t11 * t12;
t13 = t19 * t11 - t22 * t8;
t4 = t10 * t17 + t21;
t3 = t7 * t17 - t18;
t2 = t11 * t18 - t12 * t7;
t1 = t10 * t12 + t11 * t21;
t5 = [pkin(7) * t12 - t16 * t1 - t20 * t2 - t23 * t9, t13 * t12, t16 * t4 - t20 * t3, t3, 0, 0; t9 * pkin(7) + t12 * t23 + t16 * t3 + t20 * t4, t13 * t9, -t20 * t1 + t16 * t2, t1, 0, 0; 0, t22 * t11 + t15 (t16 * t10 - t20 * t7) * t8, t8 * t7, 0, 0;];
Ja_transl  = t5;
