% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPP1_jacobia_transl_2_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_transl_2_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobia_transl_2_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_transl_2_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (25->16), mult. (32->23), div. (0->0), fcn. (35->10), ass. (0->11)
t6 = sin(pkin(4));
t10 = t6 * (r_i_i_C(3) + qJ(2));
t9 = cos(qJ(1));
t8 = sin(qJ(1));
t7 = cos(pkin(6));
t5 = sin(pkin(6));
t4 = pkin(4) - pkin(6);
t3 = pkin(4) + pkin(6);
t2 = cos(t4) / 0.2e1 + cos(t3) / 0.2e1;
t1 = sin(t3) / 0.2e1 - sin(t4) / 0.2e1;
t11 = [(-t9 * t1 - t8 * t7) * r_i_i_C(1) + (-t9 * t2 + t8 * t5) * r_i_i_C(2) - t8 * pkin(1) + t9 * t10, t8 * t6, 0, 0; (-t8 * t1 + t9 * t7) * r_i_i_C(1) + (-t8 * t2 - t9 * t5) * r_i_i_C(2) + t9 * pkin(1) + t8 * t10, -t9 * t6, 0, 0; 0, cos(pkin(4)) 0, 0;];
Ja_transl  = t11;
