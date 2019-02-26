% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobia_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:32
% EndTime: 2019-02-26 22:04:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (73->15), mult. (71->18), div. (0->0), fcn. (74->6), ass. (0->17)
t32 = r_i_i_C(3) + qJ(4);
t13 = qJ(2) + qJ(3);
t10 = sin(t13);
t11 = cos(t13);
t27 = pkin(3) + r_i_i_C(1);
t19 = t32 * t10 + t27 * t11;
t31 = cos(qJ(2)) * pkin(2) + t19;
t30 = t32 * t11;
t28 = pkin(1) + t31;
t15 = sin(qJ(1));
t26 = t30 * t15;
t16 = cos(qJ(1));
t25 = t30 * t16;
t23 = r_i_i_C(2) + pkin(8) + pkin(7);
t20 = t27 * t10;
t18 = -sin(qJ(2)) * pkin(2) - t20;
t1 = [-t28 * t15 + t23 * t16, t18 * t16 + t25, -t16 * t20 + t25, t16 * t10, 0, 0; t23 * t15 + t28 * t16, t18 * t15 + t26, -t15 * t20 + t26, t15 * t10, 0, 0; 0, t31, t19, -t11, 0, 0;];
Ja_transl  = t1;
