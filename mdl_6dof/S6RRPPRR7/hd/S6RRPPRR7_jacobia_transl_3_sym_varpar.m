% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR7
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR7_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_transl_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:56
% EndTime: 2019-02-26 21:31:56
% DurationCPUTime: 0.07s
% Computational Cost: add. (40->17), mult. (89->25), div. (0->0), fcn. (110->6), ass. (0->18)
t8 = sin(qJ(2));
t9 = sin(qJ(1));
t18 = t8 * t9;
t17 = pkin(2) + r_i_i_C(1);
t11 = cos(qJ(1));
t16 = t11 * t8;
t10 = cos(qJ(2));
t15 = t9 * t10;
t14 = t10 * t11;
t13 = r_i_i_C(3) + qJ(3);
t6 = sin(pkin(6));
t12 = (pkin(8) + r_i_i_C(2)) * t6;
t7 = cos(pkin(6));
t4 = -t7 * t18 + t14;
t3 = t7 * t15 + t16;
t2 = t7 * t16 + t15;
t1 = -t7 * t14 + t18;
t5 = [-t9 * pkin(1) - t13 * t1 + t11 * t12 - t17 * t2, t13 * t4 - t17 * t3, t3, 0, 0, 0; t11 * pkin(1) + t9 * t12 + t13 * t3 + t17 * t4, -t17 * t1 + t13 * t2, t1, 0, 0, 0; 0 (t17 * t10 + t13 * t8) * t6, -t6 * t10, 0, 0, 0;];
Ja_transl  = t5;
