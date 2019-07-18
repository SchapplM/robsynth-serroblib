% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobia_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobia_transl_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:34
% EndTime: 2019-07-18 13:30:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (87->11), mult. (52->14), div. (0->0), fcn. (52->8), ass. (0->14)
t23 = pkin(6) + r_i_i_C(3);
t13 = sin(qJ(5));
t14 = cos(qJ(5));
t22 = t14 * r_i_i_C(1) - t13 * r_i_i_C(2);
t12 = qJ(2) + qJ(3);
t19 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t11 = qJ(4) + t12;
t8 = sin(t11);
t9 = cos(t11);
t18 = -t22 * t8 + t23 * t9;
t17 = t22 * t9 + t23 * t8;
t16 = t17 + pkin(3) * cos(t12);
t15 = -pkin(3) * sin(t12) + t18;
t1 = [0, -sin(qJ(2)) * pkin(2) + t15, t15, t18, t19 * t9; 0, cos(qJ(2)) * pkin(2) + t16, t16, t17, t19 * t8; 1, 0, 0, 0, t22;];
Ja_transl  = t1;
