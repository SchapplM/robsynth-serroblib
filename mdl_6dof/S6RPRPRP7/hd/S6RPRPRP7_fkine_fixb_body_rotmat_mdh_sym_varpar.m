% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:49
	% EndTime: 2020-11-04 21:37:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:49
	% EndTime: 2020-11-04 21:37:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:49
	% EndTime: 2020-11-04 21:37:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [0, -t59, t58, t59 * pkin(1) + t58 * qJ(2) + 0; 0, -t58, -t59, t58 * pkin(1) - t59 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:49
	% EndTime: 2020-11-04 21:37:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t64 = pkin(1) + pkin(7);
	t63 = cos(qJ(1));
	t62 = cos(qJ(3));
	t61 = sin(qJ(1));
	t60 = sin(qJ(3));
	t1 = [t61 * t60, t61 * t62, t63, t61 * qJ(2) + t64 * t63 + 0; -t63 * t60, -t63 * t62, t61, -t63 * qJ(2) + t64 * t61 + 0; t62, -t60, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:49
	% EndTime: 2020-11-04 21:37:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t71 = cos(qJ(1));
	t70 = sin(qJ(1));
	t69 = qJ(3) + pkin(9);
	t68 = pkin(1) + pkin(7) + qJ(4);
	t67 = cos(t69);
	t66 = sin(t69);
	t65 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t70 * t66, t70 * t67, t71, t65 * t70 + t68 * t71 + 0; -t71 * t66, -t71 * t67, t70, -t65 * t71 + t68 * t70 + 0; t67, -t66, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:50
	% EndTime: 2020-11-04 21:37:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t80 = sin(qJ(5));
	t82 = sin(qJ(1));
	t90 = t82 * t80;
	t83 = cos(qJ(5));
	t89 = t82 * t83;
	t85 = cos(qJ(1));
	t88 = t85 * t80;
	t87 = t85 * t83;
	t78 = sin(pkin(9));
	t79 = cos(pkin(9));
	t84 = cos(qJ(3));
	t86 = (t79 * pkin(4) + t78 * pkin(8) + pkin(3)) * sin(qJ(3)) - (-t78 * pkin(4) + t79 * pkin(8)) * t84 + qJ(2);
	t77 = qJ(3) + pkin(9);
	t76 = pkin(1) + pkin(7) + qJ(4);
	t75 = cos(t77);
	t74 = sin(t77);
	t1 = [t74 * t89 + t88, -t74 * t90 + t87, -t82 * t75, t76 * t85 + t86 * t82 + 0; -t74 * t87 + t90, t74 * t88 + t89, t85 * t75, t76 * t82 - t86 * t85 + 0; t75 * t83, -t75 * t80, t74, t84 * pkin(3) + t75 * pkin(4) + t74 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:50
	% EndTime: 2020-11-04 21:37:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->24), mult. (40->24), div. (0->0), fcn. (53->8), ass. (0->16)
	t98 = sin(qJ(5));
	t99 = sin(qJ(1));
	t107 = t99 * t98;
	t101 = cos(qJ(1));
	t106 = t101 * t98;
	t100 = cos(qJ(5));
	t105 = t99 * t100;
	t104 = t101 * t100;
	t103 = pkin(5) * t98 + pkin(1) + pkin(7) + qJ(4);
	t92 = t100 * pkin(5) + pkin(4);
	t96 = qJ(3) + pkin(9);
	t93 = sin(t96);
	t94 = cos(t96);
	t97 = -qJ(6) - pkin(8);
	t102 = t92 * t93 + t94 * t97 + sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t93 * t105 + t106, -t93 * t107 + t104, -t99 * t94, t103 * t101 + t102 * t99 + 0; -t93 * t104 + t107, t93 * t106 + t105, t101 * t94, -t102 * t101 + t103 * t99 + 0; t94 * t100, -t94 * t98, t93, t94 * t92 - t93 * t97 + cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end