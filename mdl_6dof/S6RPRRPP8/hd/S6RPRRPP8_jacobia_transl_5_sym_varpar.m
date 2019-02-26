% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobia_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:20
% EndTime: 2019-02-26 21:00:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (64->26), mult. (144->33), div. (0->0), fcn. (167->6), ass. (0->23)
t11 = cos(qJ(3));
t19 = pkin(8) + r_i_i_C(1);
t10 = cos(qJ(4));
t15 = r_i_i_C(3) + qJ(5);
t20 = pkin(4) - r_i_i_C(2);
t7 = sin(qJ(4));
t23 = t20 * t10 + t15 * t7 + pkin(3);
t8 = sin(qJ(3));
t26 = t23 * t11 + t19 * t8;
t24 = t19 * t11;
t9 = sin(qJ(1));
t22 = t9 * t7;
t21 = pkin(1) + pkin(7);
t12 = cos(qJ(1));
t18 = t12 * t7;
t17 = t9 * t10;
t16 = t12 * t10;
t14 = t8 * pkin(3) + qJ(2) - t24;
t4 = t8 * t16 - t22;
t3 = t8 * t18 + t17;
t2 = t8 * t17 + t18;
t1 = t8 * t22 - t16;
t5 = [t14 * t12 + t15 * t3 + t20 * t4 - t21 * t9, t9, t26 * t9, -t20 * t1 + t15 * t2, t1, 0; t15 * t1 + t21 * t12 + t14 * t9 + t20 * t2, -t12, -t26 * t12, -t15 * t4 + t20 * t3, -t3, 0; 0, 0, -t23 * t8 + t24 (t15 * t10 - t20 * t7) * t11, t11 * t7, 0;];
Ja_transl  = t5;
