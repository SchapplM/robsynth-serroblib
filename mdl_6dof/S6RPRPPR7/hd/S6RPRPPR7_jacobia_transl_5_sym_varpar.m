% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:33
% EndTime: 2019-02-26 20:42:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (55->14), mult. (55->14), div. (0->0), fcn. (62->6), ass. (0->12)
t3 = qJ(3) + pkin(9);
t1 = sin(t3);
t11 = r_i_i_C(3) + qJ(5);
t12 = pkin(4) - r_i_i_C(2);
t2 = cos(t3);
t16 = -t12 * t1 + t11 * t2 - sin(qJ(3)) * pkin(3);
t15 = t11 * t1 + t12 * t2 + cos(qJ(3)) * pkin(3);
t10 = pkin(1) + r_i_i_C(1) + qJ(4) + pkin(7);
t9 = qJ(2) - t16;
t8 = cos(qJ(1));
t6 = sin(qJ(1));
t4 = [-t10 * t6 + t9 * t8, t6, t15 * t6, t8, -t6 * t2, 0; t10 * t8 + t9 * t6, -t8, -t15 * t8, t6, t8 * t2, 0; 0, 0, t16, 0, t1, 0;];
Ja_transl  = t4;
