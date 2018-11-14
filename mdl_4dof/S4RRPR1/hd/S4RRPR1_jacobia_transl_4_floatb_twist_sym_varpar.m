% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRPR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S4RRPR1_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR1_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:40
% EndTime: 2018-11-14 13:53:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (64->11), mult. (22->10), div. (0->0), fcn. (22->8), ass. (0->10)
t10 = qJ(1) + qJ(2);
t8 = pkin(7) + t10;
t7 = qJ(4) + t8;
t3 = sin(t7);
t4 = cos(t7);
t14 = r_i_i_C(1) * t4 - r_i_i_C(2) * t3;
t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
t12 = t14 + pkin(3) * cos(t8) + pkin(2) * cos(t10);
t11 = -pkin(2) * sin(t10) - pkin(3) * sin(t8) + t13;
t1 = [-sin(qJ(1)) * pkin(1) + t11, t11, 0, t13; cos(qJ(1)) * pkin(1) + t12, t12, 0, t14; 0, 0, 1, 0;];
Ja_transl  = t1;
