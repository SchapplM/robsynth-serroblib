% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:53
% EndTime: 2019-02-26 20:56:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (158->27), mult. (181->39), div. (0->0), fcn. (209->8), ass. (0->21)
t13 = cos(qJ(3));
t11 = sin(qJ(3));
t17 = pkin(8) - r_i_i_C(3) - qJ(6);
t15 = t17 * t11;
t23 = t13 * pkin(3) + pkin(2) + t15;
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t18 = pkin(4) + pkin(5) + r_i_i_C(1);
t19 = r_i_i_C(2) + qJ(5);
t22 = t19 * t10 + t18 * t12 + pkin(3);
t21 = t10 * t13;
t20 = t12 * t13;
t14 = -t22 * t11 + t17 * t13;
t9 = qJ(1) + pkin(9);
t8 = cos(t9);
t7 = sin(t9);
t4 = t7 * t10 + t8 * t20;
t3 = -t7 * t12 + t8 * t21;
t2 = -t8 * t10 + t7 * t20;
t1 = t8 * t12 + t7 * t21;
t5 = [-sin(qJ(1)) * pkin(1) + t8 * pkin(7) - t19 * t1 - t18 * t2 - t23 * t7, 0, t14 * t8, -t18 * t3 + t19 * t4, t3, -t8 * t11; cos(qJ(1)) * pkin(1) + t7 * pkin(7) + t19 * t3 + t18 * t4 + t23 * t8, 0, t14 * t7, -t18 * t1 + t19 * t2, t1, -t7 * t11; 0, 1, t22 * t13 + t15 (-t18 * t10 + t19 * t12) * t11, t11 * t10, t13;];
Ja_transl  = t5;
