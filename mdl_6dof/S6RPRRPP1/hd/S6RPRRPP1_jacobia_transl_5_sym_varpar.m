% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:16
% EndTime: 2019-02-26 20:56:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (123->32), mult. (108->47), div. (0->0), fcn. (119->10), ass. (0->26)
t16 = cos(qJ(3));
t14 = sin(qJ(3));
t23 = r_i_i_C(3) + qJ(5) + pkin(8);
t19 = t23 * t14;
t15 = cos(qJ(4));
t5 = pkin(4) * t15 + pkin(3);
t27 = t16 * t5 + pkin(2) + t19;
t13 = sin(qJ(4));
t26 = pkin(4) * t13;
t11 = qJ(1) + pkin(9);
t7 = sin(t11);
t25 = t16 * t7;
t9 = cos(t11);
t24 = t16 * t9;
t22 = t13 * t16;
t20 = pkin(7) + t26;
t10 = qJ(4) + pkin(10);
t6 = sin(t10);
t8 = cos(t10);
t18 = r_i_i_C(1) * t8 - r_i_i_C(2) * t6 + t5;
t17 = -t18 * t14 + t23 * t16;
t4 = t8 * t24 + t7 * t6;
t3 = -t6 * t24 + t7 * t8;
t2 = -t8 * t25 + t9 * t6;
t1 = t6 * t25 + t9 * t8;
t12 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t20 * t9 - t27 * t7, 0, t17 * t9, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (t15 * t7 - t9 * t22) * pkin(4), t9 * t14, 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t20 * t7 + t27 * t9, 0, t17 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t15 * t9 - t7 * t22) * pkin(4), t7 * t14, 0; 0, 1, t18 * t16 + t19 (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8 - t26) * t14, -t16, 0;];
Ja_transl  = t12;
