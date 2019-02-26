% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:27:33
% EndTime: 2019-02-26 20:27:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->18), mult. (73->20), div. (0->0), fcn. (95->8), ass. (0->15)
t23 = pkin(1) + pkin(2);
t9 = qJ(4) + pkin(10);
t7 = sin(t9);
t8 = cos(t9);
t22 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2) - cos(qJ(4)) * pkin(4);
t20 = r_i_i_C(3) + qJ(5) + pkin(7);
t19 = cos(qJ(1));
t18 = sin(qJ(1));
t17 = cos(pkin(9));
t16 = sin(pkin(9));
t14 = pkin(3) - t22;
t13 = sin(qJ(4)) * pkin(4) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
t2 = t19 * t16 - t18 * t17;
t1 = -t18 * t16 - t19 * t17;
t3 = [t19 * qJ(2) + t20 * t1 + t14 * t2 - t23 * t18, t18, 0, t13 * t1, t2, 0; t18 * qJ(2) - t14 * t1 + t23 * t19 + t20 * t2, -t19, 0, t13 * t2, -t1, 0; 0, 0, -1, t22, 0, 0;];
Ja_transl  = t3;
