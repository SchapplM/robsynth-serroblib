% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP4_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP4_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobia_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:46
% EndTime: 2019-02-26 21:26:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (45->20), mult. (104->30), div. (0->0), fcn. (120->6), ass. (0->18)
t15 = r_i_i_C(2) + qJ(3);
t7 = sin(qJ(2));
t12 = t15 * t7;
t9 = cos(qJ(2));
t20 = t9 * pkin(2) + pkin(1) + t12;
t14 = r_i_i_C(3) + qJ(4);
t17 = pkin(3) + r_i_i_C(1);
t5 = sin(pkin(9));
t6 = cos(pkin(9));
t19 = t14 * t5 + t17 * t6 + pkin(2);
t8 = sin(qJ(1));
t18 = t8 * t9;
t10 = cos(qJ(1));
t16 = t10 * t9;
t11 = t15 * t9 - t19 * t7;
t3 = t5 * t16 - t8 * t6;
t1 = t10 * t6 + t5 * t18;
t2 = [pkin(7) * t10 + t17 * (t10 * t5 - t6 * t18) - t14 * t1 - t20 * t8, t11 * t10, t10 * t7, t3, 0, 0; t8 * pkin(7) + t17 * (t6 * t16 + t8 * t5) + t14 * t3 + t20 * t10, t11 * t8, t8 * t7, t1, 0, 0; 0, t19 * t9 + t12, -t9, t7 * t5, 0, 0;];
Ja_transl  = t2;
