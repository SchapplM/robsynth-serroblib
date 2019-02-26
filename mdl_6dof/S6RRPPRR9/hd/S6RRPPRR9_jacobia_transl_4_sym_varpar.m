% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR9_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->18), mult. (118->28), div. (0->0), fcn. (148->6), ass. (0->16)
t17 = r_i_i_C(2) + qJ(3);
t16 = cos(pkin(6));
t15 = pkin(2) + r_i_i_C(3) + qJ(4);
t9 = sin(qJ(1));
t14 = t9 * t16;
t11 = cos(qJ(1));
t13 = t11 * t16;
t7 = sin(pkin(6));
t12 = (pkin(3) + pkin(8) + r_i_i_C(1)) * t7;
t10 = cos(qJ(2));
t8 = sin(qJ(2));
t4 = t10 * t11 - t8 * t14;
t3 = t10 * t14 + t11 * t8;
t2 = t9 * t10 + t8 * t13;
t1 = -t10 * t13 + t8 * t9;
t5 = [-t9 * pkin(1) - t17 * t1 + t11 * t12 - t15 * t2, -t15 * t3 + t17 * t4, t3, t4, 0, 0; pkin(1) * t11 + t9 * t12 + t15 * t4 + t17 * t3, -t15 * t1 + t17 * t2, t1, t2, 0, 0; 0 (t15 * t10 + t17 * t8) * t7, -t7 * t10, t7 * t8, 0, 0;];
Ja_transl  = t5;
