% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPRP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPRP1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_jacobia_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPRP1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_jacobia_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:32:12
% EndTime: 2019-02-26 19:32:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (55->11), mult. (20->8), div. (0->0), fcn. (22->6), ass. (0->9)
t14 = pkin(3) + r_i_i_C(1);
t13 = qJ(4) + r_i_i_C(3);
t10 = qJ(1) + pkin(6);
t9 = qJ(3) + t10;
t7 = sin(t9);
t8 = cos(t9);
t12 = t13 * t7 + t14 * t8;
t11 = t13 * t8 - t14 * t7;
t1 = [-pkin(2) * sin(t10) - sin(qJ(1)) * pkin(1) + t11, 0, t11, t7; pkin(2) * cos(t10) + cos(qJ(1)) * pkin(1) + t12, 0, t12, -t8; 0, 1, 0, 0;];
Ja_transl  = t1;
