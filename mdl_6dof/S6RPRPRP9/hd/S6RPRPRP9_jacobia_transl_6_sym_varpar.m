% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:11
% EndTime: 2019-02-26 20:48:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (130->30), mult. (155->39), div. (0->0), fcn. (181->8), ass. (0->23)
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t22 = r_i_i_C(2) + pkin(8) + qJ(4);
t19 = r_i_i_C(3) + qJ(6);
t23 = pkin(5) + r_i_i_C(1);
t6 = cos(pkin(9)) * pkin(4) + pkin(3);
t9 = pkin(9) + qJ(5);
t7 = sin(t9);
t8 = cos(t9);
t24 = t19 * t7 + t23 * t8 + t6;
t27 = t22 * t12 + t24 * t14;
t25 = t22 * t14;
t13 = sin(qJ(1));
t21 = t12 * t13;
t15 = cos(qJ(1));
t20 = t12 * t15;
t18 = pkin(4) * sin(pkin(9)) + pkin(1) + pkin(7);
t17 = t12 * t6 + qJ(2) - t25;
t4 = -t13 * t7 + t8 * t20;
t3 = t13 * t8 + t7 * t20;
t2 = t15 * t7 + t8 * t21;
t1 = -t15 * t8 + t7 * t21;
t5 = [-t18 * t13 + t17 * t15 + t19 * t3 + t23 * t4, t13, t27 * t13, -t13 * t14, -t23 * t1 + t19 * t2, t1; t19 * t1 + t17 * t13 + t18 * t15 + t23 * t2, -t15, -t27 * t15, t15 * t14, -t19 * t4 + t23 * t3, -t3; 0, 0, -t12 * t24 + t25, t12 (t19 * t8 - t23 * t7) * t14, t14 * t7;];
Ja_transl  = t5;
