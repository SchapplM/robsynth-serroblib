% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobia_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobia_transl_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:50
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (69->10), mult. (55->16), div. (0->0), fcn. (55->8), ass. (0->13)
t11 = qJ(3) + qJ(4);
t7 = sin(t11);
t9 = cos(t11);
t19 = t9 * r_i_i_C(1) - t7 * r_i_i_C(2);
t23 = t19 + cos(qJ(3)) * pkin(2);
t18 = -r_i_i_C(1) * t7 - r_i_i_C(2) * t9;
t12 = qJ(1) + qJ(2);
t10 = cos(t12);
t8 = sin(t12);
t17 = t8 * r_i_i_C(3) + t10 * t23;
t16 = -sin(qJ(3)) * pkin(2) + t18;
t15 = t10 * r_i_i_C(3) - t23 * t8;
t1 = [-sin(qJ(1)) * pkin(1) + t15, t15, t16 * t10, t18 * t10, 0; cos(qJ(1)) * pkin(1) + t17, t17, t16 * t8, t18 * t8, 0; 0, 0, t23, t19, 0;];
Ja_transl  = t1;
