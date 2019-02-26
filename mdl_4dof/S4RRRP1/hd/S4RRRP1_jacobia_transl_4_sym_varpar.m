% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RRRP1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_jacobia_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRRP1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_jacobia_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:37:00
% EndTime: 2019-02-26 19:37:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (58->9), mult. (24->8), div. (0->0), fcn. (24->6), ass. (0->10)
t13 = pkin(3) + r_i_i_C(1);
t8 = qJ(1) + qJ(2);
t7 = qJ(3) + t8;
t4 = sin(t7);
t5 = cos(t7);
t12 = -t4 * r_i_i_C(2) + t13 * t5;
t11 = t12 + pkin(2) * cos(t8);
t10 = -t5 * r_i_i_C(2) - t13 * t4;
t9 = -pkin(2) * sin(t8) + t10;
t1 = [-sin(qJ(1)) * pkin(1) + t9, t9, t10, 0; cos(qJ(1)) * pkin(1) + t11, t11, t12, 0; 0, 0, 0, 1;];
Ja_transl  = t1;
