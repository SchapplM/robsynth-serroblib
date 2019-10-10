% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S2RR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% JaD_transl [3x2]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S2RR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_jacobiaD_transl_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_jacobiaD_transl_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S2RR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S2RR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_jacobiaD_transl_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0; 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t66 = sin(qJ(2));
	t68 = cos(qJ(2));
	t70 = qJD(2) * (r_i_i_C(1) * t66 + r_i_i_C(2) * t68);
	t77 = -pkin(1) - r_i_i_C(3);
	t67 = sin(qJ(1));
	t76 = qJD(1) * t67;
	t69 = cos(qJ(1));
	t75 = qJD(1) * t69;
	t74 = qJD(2) * t67;
	t73 = qJD(2) * t69;
	t72 = r_i_i_C(1) * t68 - r_i_i_C(2) * t66;
	t1 = [t67 * t70 + (t77 * t67 - t72 * t69) * qJD(1), (t66 * t73 + t68 * t76) * r_i_i_C(2) + (t66 * t76 - t68 * t73) * r_i_i_C(1); 0, -t70; t69 * t70 + (t72 * t67 + t77 * t69) * qJD(1), (-t66 * t74 + t68 * t75) * r_i_i_C(2) + (t66 * t75 + t68 * t74) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,2);
end