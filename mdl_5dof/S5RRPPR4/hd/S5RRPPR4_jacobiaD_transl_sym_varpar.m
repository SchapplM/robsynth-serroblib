% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (58->10), mult. (42->12), div. (0->0), fcn. (24->4), ass. (0->12)
	t31 = qJ(3) + r_i_i_C(3);
	t30 = -pkin(2) - r_i_i_C(1);
	t24 = qJ(1) + qJ(2);
	t21 = sin(t24);
	t23 = qJD(1) + qJD(2);
	t29 = t23 * t21;
	t22 = cos(t24);
	t28 = t23 * t22;
	t27 = pkin(1) * qJD(1);
	t26 = t21 * qJD(3) + t28 * t31 + t30 * t29;
	t25 = t22 * qJD(3) + (-t21 * t31 + t30 * t22) * t23;
	t1 = [-cos(qJ(1)) * t27 + t25, t25, t28, 0, 0; -sin(qJ(1)) * t27 + t26, t26, t29, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (94->15), mult. (74->22), div. (0->0), fcn. (52->6), ass. (0->15)
	t90 = -pkin(2) - pkin(3);
	t82 = qJ(1) + qJ(2);
	t79 = sin(t82);
	t81 = qJD(1) + qJD(2);
	t89 = t81 * t79;
	t80 = cos(t82);
	t88 = t81 * t80;
	t87 = pkin(1) * qJD(1);
	t83 = sin(pkin(8));
	t84 = cos(pkin(8));
	t74 = (t79 * t83 + t80 * t84) * t81;
	t75 = (-t79 * t84 + t80 * t83) * t81;
	t86 = t75 * r_i_i_C(1) + t74 * r_i_i_C(2) + qJ(3) * t88 + t79 * qJD(3) + t90 * t89;
	t85 = -t74 * r_i_i_C(1) + t75 * r_i_i_C(2) + t80 * qJD(3) + (-qJ(3) * t79 + t90 * t80) * t81;
	t1 = [-cos(qJ(1)) * t87 + t85, t85, t88, 0, 0; -sin(qJ(1)) * t87 + t86, t86, t89, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:05
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (237->27), mult. (230->43), div. (0->0), fcn. (190->8), ass. (0->26)
	t97 = -pkin(7) - r_i_i_C(3);
	t76 = qJ(1) + qJ(2);
	t73 = sin(t76);
	t74 = cos(t76);
	t77 = sin(pkin(8));
	t78 = cos(pkin(8));
	t67 = -t73 * t78 + t74 * t77;
	t75 = qJD(1) + qJD(2);
	t66 = t67 * t75;
	t68 = t73 * t77 + t74 * t78;
	t80 = cos(qJ(5));
	t79 = sin(qJ(5));
	t88 = qJD(5) * t79;
	t96 = -t66 * t80 + t68 * t88;
	t65 = t68 * t75;
	t95 = -t65 * t80 - t67 * t88;
	t94 = -pkin(2) - pkin(3);
	t91 = t75 * t73;
	t90 = t75 * t74;
	t89 = pkin(1) * qJD(1);
	t87 = qJD(5) * t80;
	t84 = -t65 * t79 + t67 * t87;
	t83 = -t66 * t79 - t68 * t87;
	t82 = t66 * pkin(4) - t96 * r_i_i_C(1) + t83 * r_i_i_C(2) + qJ(3) * t90 + t73 * qJD(3) + t97 * t65 + t94 * t91;
	t81 = -t65 * pkin(4) + t74 * qJD(3) - t84 * r_i_i_C(2) + (-qJ(3) * t73 + t94 * t74) * t75 + t97 * t66 + t95 * r_i_i_C(1);
	t1 = [-cos(qJ(1)) * t89 + t81, t81, t90, 0, t83 * r_i_i_C(1) + t96 * r_i_i_C(2); -sin(qJ(1)) * t89 + t82, t82, t91, 0, t84 * r_i_i_C(1) + t95 * r_i_i_C(2); 0, 0, 0, 0, (r_i_i_C(1) * t79 + r_i_i_C(2) * t80) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end