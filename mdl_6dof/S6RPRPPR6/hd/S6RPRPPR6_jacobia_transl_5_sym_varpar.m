% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:58
% EndTime: 2019-02-26 20:41:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (67->18), mult. (76->18), div. (0->0), fcn. (87->8), ass. (0->14)
t3 = qJ(3) + pkin(9);
t1 = sin(t3);
t4 = sin(pkin(10));
t5 = cos(pkin(10));
t13 = r_i_i_C(1) * t5 - r_i_i_C(2) * t4 + pkin(4);
t14 = r_i_i_C(3) + qJ(5);
t2 = cos(t3);
t18 = -t13 * t1 + t14 * t2 - sin(qJ(3)) * pkin(3);
t17 = t14 * t1 + t13 * t2 + cos(qJ(3)) * pkin(3);
t12 = t4 * r_i_i_C(1) + r_i_i_C(2) * t5 + pkin(1) + pkin(7) + qJ(4);
t11 = qJ(2) - t18;
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t6 = [t11 * t10 - t12 * t8, t8, t17 * t8, t10, -t8 * t2, 0; t12 * t10 + t11 * t8, -t10, -t17 * t10, t8, t10 * t2, 0; 0, 0, t18, 0, t1, 0;];
Ja_transl  = t6;
