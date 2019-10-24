% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRPR2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->3), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->5)
	t14 = pkin(8) + qJ(3);
	t12 = sin(t14);
	t13 = cos(t14);
	t17 = qJD(3) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, 0, cos(pkin(7)) * t17, 0, 0; 0, 0, sin(pkin(7)) * t17, 0, 0; 0, 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->7), mult. (42->14), div. (0->0), fcn. (30->4), ass. (0->10)
	t102 = -pkin(3) + r_i_i_C(2);
	t101 = r_i_i_C(3) + qJ(4);
	t96 = pkin(8) + qJ(3);
	t95 = cos(t96);
	t100 = qJD(3) * t95;
	t94 = sin(t96);
	t99 = qJD(4) * t95 + (-t101 * t94 + t102 * t95) * qJD(3);
	t98 = cos(pkin(7));
	t97 = sin(pkin(7));
	t1 = [0, 0, t99 * t98, t98 * t100, 0; 0, 0, t99 * t97, t97 * t100, 0; 0, 0, t94 * qJD(4) + (t101 * t95 + t102 * t94) * qJD(3), qJD(3) * t94, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (79->21), mult. (134->43), div. (0->0), fcn. (102->6), ass. (0->21)
	t138 = sin(pkin(7));
	t140 = sin(qJ(5));
	t154 = t138 * t140;
	t141 = cos(qJ(5));
	t153 = t138 * t141;
	t139 = cos(pkin(7));
	t152 = t139 * t140;
	t151 = t139 * t141;
	t137 = pkin(8) + qJ(3);
	t135 = sin(t137);
	t150 = qJD(3) * t135;
	t136 = cos(t137);
	t149 = qJD(3) * t136;
	t148 = qJD(5) * t136;
	t147 = -pkin(3) - pkin(6) - r_i_i_C(3);
	t146 = r_i_i_C(1) * t141 - r_i_i_C(2) * t140;
	t145 = r_i_i_C(1) * t140 + r_i_i_C(2) * t141 + qJ(4);
	t144 = t146 * t149;
	t143 = t146 * qJD(5) + qJD(4);
	t142 = t143 * t136 + (-t145 * t135 + t147 * t136) * qJD(3);
	t1 = [0, 0, t142 * t139, t139 * t149, t139 * t144 + ((-t135 * t152 - t153) * r_i_i_C(1) + (-t135 * t151 + t154) * r_i_i_C(2)) * qJD(5); 0, 0, t142 * t138, t138 * t149, t138 * t144 + ((-t135 * t154 + t151) * r_i_i_C(1) + (-t135 * t153 - t152) * r_i_i_C(2)) * qJD(5); 0, 0, t143 * t135 + (t147 * t135 + t145 * t136) * qJD(3), t150, (-t140 * t150 + t141 * t148) * r_i_i_C(2) + (t140 * t148 + t141 * t150) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end