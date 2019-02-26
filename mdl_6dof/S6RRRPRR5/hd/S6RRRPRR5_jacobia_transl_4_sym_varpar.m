% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR5_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:22
% EndTime: 2019-02-26 22:18:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (73->17), mult. (71->20), div. (0->0), fcn. (74->6), ass. (0->17)
t33 = r_i_i_C(3) + qJ(4);
t14 = qJ(2) + qJ(3);
t11 = sin(t14);
t12 = cos(t14);
t19 = (pkin(3) - r_i_i_C(2)) * t12 + t33 * t11;
t32 = cos(qJ(2)) * pkin(2) + t19;
t31 = t33 * t12;
t30 = pkin(1) + t32;
t27 = r_i_i_C(1) + pkin(8) + pkin(7);
t16 = sin(qJ(1));
t26 = t16 * t11;
t17 = cos(qJ(1));
t25 = t17 * t11;
t22 = r_i_i_C(2) * t26 + t31 * t16;
t21 = r_i_i_C(2) * t25 + t31 * t17;
t20 = -sin(qJ(2)) * pkin(2) - pkin(3) * t11;
t1 = [-t30 * t16 + t27 * t17, t20 * t17 + t21, -pkin(3) * t25 + t21, t25, 0, 0; t27 * t16 + t30 * t17, t20 * t16 + t22, -pkin(3) * t26 + t22, t26, 0, 0; 0, t32, t19, -t12, 0, 0;];
Ja_transl  = t1;
