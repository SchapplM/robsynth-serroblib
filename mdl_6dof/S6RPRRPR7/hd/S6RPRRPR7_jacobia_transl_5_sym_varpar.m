% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:04:21
% EndTime: 2019-02-26 21:04:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (74->19), mult. (55->18), div. (0->0), fcn. (59->8), ass. (0->17)
t11 = qJ(3) + qJ(4);
t7 = pkin(10) + t11;
t5 = sin(t7);
t6 = cos(t7);
t15 = -t5 * r_i_i_C(1) - t6 * r_i_i_C(2) - pkin(4) * sin(t11);
t24 = t15 - sin(qJ(3)) * pkin(3);
t22 = pkin(4) * cos(t11);
t21 = r_i_i_C(1) * t6;
t20 = r_i_i_C(2) * t5;
t18 = pkin(1) + r_i_i_C(3) + qJ(5) + pkin(8) + pkin(7);
t16 = qJ(2) - t24;
t14 = cos(qJ(1));
t13 = sin(qJ(1));
t4 = t14 * t20;
t3 = t13 * t21;
t2 = t22 + cos(qJ(3)) * pkin(3);
t1 = [-t18 * t13 + t16 * t14, t13, t3 + (t2 - t20) * t13, t3 + (-t20 + t22) * t13, t14, 0; t16 * t13 + t18 * t14, -t14, t4 + (-t2 - t21) * t14, t4 + (-t21 - t22) * t14, t13, 0; 0, 0, t24, t15, 0, 0;];
Ja_transl  = t1;
