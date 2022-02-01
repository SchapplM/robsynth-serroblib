% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(8);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t42 = qJ(1) + pkin(8);
	t40 = qJ(3) + t42;
	t38 = sin(t40);
	t39 = cos(t40);
	t41 = qJD(1) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t41;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t41;
	t1 = [(-cos(t42) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t44, 0, t44, 0, 0; t43 + (-sin(t42) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t43, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (100->14), mult. (62->17), div. (0->0), fcn. (38->8), ass. (0->13)
	t39 = qJ(1) + pkin(8);
	t37 = qJ(3) + t39;
	t35 = sin(t37);
	t38 = qJD(1) + qJD(3);
	t47 = t38 * t35;
	t49 = qJ(4) + r_i_i_C(3);
	t48 = qJD(4) + r_i_i_C(2) * t38 * sin(pkin(9));
	t36 = cos(t37);
	t46 = t38 * t36;
	t44 = -r_i_i_C(1) * cos(pkin(9)) - pkin(3);
	t43 = t48 * t35 + t44 * t47 + t49 * t46;
	t42 = -t49 * t47 + (t44 * t38 + t48) * t36;
	t1 = [(-cos(t39) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t42, 0, t42, t46, 0; (-sin(t39) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t43, 0, t43, t47, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (267->29), mult. (204->47), div. (0->0), fcn. (166->10), ass. (0->25)
	t175 = qJ(1) + pkin(8);
	t173 = qJ(3) + t175;
	t171 = sin(t173);
	t174 = qJD(1) + qJD(3);
	t190 = t174 * t171;
	t172 = cos(t173);
	t189 = t174 * t172;
	t177 = cos(pkin(9));
	t178 = sin(qJ(5));
	t188 = t177 * t178;
	t179 = cos(qJ(5));
	t187 = t177 * t179;
	t186 = -t171 * t178 - t172 * t187;
	t185 = -t171 * t179 + t172 * t188;
	t184 = t171 * t187 - t172 * t178;
	t183 = t171 * t188 + t172 * t179;
	t176 = sin(pkin(9));
	t182 = -pkin(4) * t177 - pkin(3) + (-pkin(7) - r_i_i_C(3)) * t176;
	t164 = t186 * qJD(5) + t183 * t174;
	t165 = t185 * qJD(5) + t184 * t174;
	t181 = -t165 * r_i_i_C(1) + t164 * r_i_i_C(2) + qJ(4) * t189 + t171 * qJD(4) + t182 * t190;
	t166 = t184 * qJD(5) + t185 * t174;
	t167 = t183 * qJD(5) + t186 * t174;
	t180 = t166 * r_i_i_C(2) + t167 * r_i_i_C(1) + t172 * qJD(4) + (-qJ(4) * t171 + t182 * t172) * t174;
	t1 = [(-cos(t175) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t180, 0, t180, t189, t164 * r_i_i_C(1) + t165 * r_i_i_C(2); (-sin(t175) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t181, 0, t181, t190, -t166 * r_i_i_C(1) + t167 * r_i_i_C(2); 0, 0, 0, 0, (-r_i_i_C(1) * t179 + r_i_i_C(2) * t178) * t176 * qJD(5);];
	JaD_transl = t1;
end