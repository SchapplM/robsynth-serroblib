% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S4PRRR1_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR1_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:45
% EndTime: 2018-11-14 13:44:45
% DurationCPUTime: 0.07s
% Computational Cost: add. (58->9), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->10)
t8 = pkin(7) + qJ(2);
t7 = qJ(3) + t8;
t6 = qJ(4) + t7;
t3 = sin(t6);
t4 = cos(t6);
t12 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
t11 = t12 + pkin(3) * cos(t7);
t10 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
t9 = -pkin(3) * sin(t7) + t10;
t1 = [0, -pkin(2) * sin(t8) + t9, t9, t10; 0, pkin(2) * cos(t8) + t11, t11, t12; 1, 0, 0, 0;];
Ja_transl  = t1;
