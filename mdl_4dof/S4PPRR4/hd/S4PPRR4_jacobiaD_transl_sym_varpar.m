% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PPRR4
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 11:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:34
	% EndTime: 2019-12-29 11:56:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:34
	% EndTime: 2019-12-29 11:56:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:29
	% EndTime: 2019-12-29 11:56:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:34
	% EndTime: 2019-12-29 11:56:34
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->3), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->5)
	t14 = pkin(7) + qJ(3);
	t12 = sin(t14);
	t13 = cos(t14);
	t17 = qJD(3) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, 0, cos(pkin(6)) * t17, 0; 0, 0, sin(pkin(6)) * t17, 0; 0, 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(3), 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:56:35
	% EndTime: 2019-12-29 11:56:35
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (64->20), mult. (116->41), div. (0->0), fcn. (87->6), ass. (0->20)
	t162 = pkin(5) + r_i_i_C(3);
	t147 = sin(pkin(6));
	t149 = sin(qJ(4));
	t161 = t147 * t149;
	t150 = cos(qJ(4));
	t160 = t147 * t150;
	t148 = cos(pkin(6));
	t159 = t148 * t149;
	t158 = t148 * t150;
	t146 = pkin(7) + qJ(3);
	t145 = cos(t146);
	t157 = qJD(3) * t145;
	t144 = sin(t146);
	t156 = qJD(4) * t144;
	t155 = r_i_i_C(1) * t149 + r_i_i_C(2) * t150;
	t154 = -r_i_i_C(1) * t150 + r_i_i_C(2) * t149 - pkin(3);
	t153 = t155 * t144;
	t152 = qJD(3) * t153;
	t151 = qJD(4) * t153 + (-t162 * t144 + t154 * t145) * qJD(3);
	t1 = [0, 0, t151 * t148, t148 * t152 + ((-t145 * t158 - t161) * r_i_i_C(1) + (t145 * t159 - t160) * r_i_i_C(2)) * qJD(4); 0, 0, t151 * t147, t147 * t152 + ((-t145 * t160 + t159) * r_i_i_C(1) + (t145 * t161 + t158) * r_i_i_C(2)) * qJD(4); 0, 0, -t155 * t145 * qJD(4) + (t154 * t144 + t162 * t145) * qJD(3), (t149 * t156 - t150 * t157) * r_i_i_C(2) + (-t149 * t157 - t150 * t156) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end