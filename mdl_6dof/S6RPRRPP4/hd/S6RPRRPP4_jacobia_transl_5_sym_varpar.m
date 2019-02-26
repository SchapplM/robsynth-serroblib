% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP4
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
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:02
% EndTime: 2019-02-26 20:58:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (119->32), mult. (108->45), div. (0->0), fcn. (121->9), ass. (0->27)
t26 = r_i_i_C(3) + qJ(5) + pkin(8);
t11 = pkin(9) + qJ(3);
t7 = sin(t11);
t22 = t26 * t7;
t17 = cos(qJ(4));
t6 = pkin(4) * t17 + pkin(3);
t9 = cos(t11);
t30 = t22 + t9 * t6 + cos(pkin(9)) * pkin(2) + pkin(1);
t15 = sin(qJ(4));
t29 = pkin(4) * t15;
t16 = sin(qJ(1));
t28 = t16 * t9;
t18 = cos(qJ(1));
t27 = t18 * t9;
t12 = qJ(4) + pkin(10);
t10 = cos(t12);
t25 = t10 * t18;
t24 = t16 * t10;
t21 = pkin(7) + qJ(2) + t29;
t8 = sin(t12);
t20 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + t6;
t19 = -t20 * t7 + t26 * t9;
t4 = t16 * t8 + t9 * t25;
t3 = -t8 * t27 + t24;
t2 = t18 * t8 - t9 * t24;
t1 = t8 * t28 + t25;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t30 * t16 + t21 * t18, t16, t19 * t18, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t15 * t27 + t16 * t17) * pkin(4), t18 * t7, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t21 * t16 + t30 * t18, -t18, t19 * t16, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t15 * t28 - t17 * t18) * pkin(4), t16 * t7, 0; 0, 0, t20 * t9 + t22 (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10 - t29) * t7, -t9, 0;];
Ja_transl  = t5;
