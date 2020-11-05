% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t59 = cos(pkin(9));
	t58 = sin(pkin(9));
	t1 = [t61 * t59, -t61 * t58, t60, t61 * pkin(1) + t60 * qJ(2) + 0; t60 * t59, -t60 * t58, -t61, t60 * pkin(1) - t61 * qJ(2) + 0; t58, t59, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = pkin(7) + qJ(2);
	t65 = pkin(9) + qJ(3);
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t68 * t64, -t68 * t63, t67, t68 * t62 + t66 * t67 + 0; t67 * t64, -t67 * t63, -t68, t67 * t62 - t68 * t66 + 0; t63, t64, 0, sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t72 = pkin(9) + qJ(3);
	t70 = sin(t72);
	t71 = cos(t72);
	t76 = pkin(3) * t71 + qJ(4) * t70 + cos(pkin(9)) * pkin(2) + pkin(1);
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t73 = pkin(7) + qJ(2);
	t1 = [t74, -t75 * t71, t75 * t70, t73 * t74 + t76 * t75 + 0; -t75, -t74 * t71, t74 * t70, -t75 * t73 + t76 * t74 + 0; 0, -t70, -t71, t70 * pkin(3) - t71 * qJ(4) + sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (37->26), div. (0->0), fcn. (50->10), ass. (0->17)
	t82 = sin(pkin(10));
	t87 = sin(qJ(1));
	t92 = t87 * t82;
	t84 = cos(pkin(10));
	t91 = t87 * t84;
	t88 = cos(qJ(1));
	t90 = t88 * t82;
	t89 = t88 * t84;
	t86 = pkin(3) + qJ(5);
	t85 = cos(pkin(9));
	t83 = sin(pkin(9));
	t81 = pkin(9) + qJ(3);
	t80 = pkin(4) + pkin(7) + qJ(2);
	t79 = cos(t81);
	t78 = sin(t81);
	t77 = (qJ(4) * t83 + t86 * t85) * cos(qJ(3)) + (qJ(4) * t85 - t83 * t86) * sin(qJ(3)) + t85 * pkin(2) + pkin(1);
	t1 = [t78 * t90 + t91, t78 * t89 - t92, t88 * t79, t77 * t88 + t80 * t87 + 0; t78 * t92 - t89, t78 * t91 + t90, t87 * t79, t77 * t87 - t80 * t88 + 0; -t79 * t82, -t79 * t84, t78, t83 * pkin(2) - t79 * qJ(4) + t86 * t78 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:26
	% EndTime: 2020-11-04 21:34:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (71->30), mult. (51->31), div. (0->0), fcn. (58->14), ass. (0->19)
	t105 = sin(qJ(1));
	t101 = pkin(10) + qJ(6);
	t96 = sin(t101);
	t110 = t105 * t96;
	t98 = cos(t101);
	t109 = t105 * t98;
	t106 = cos(qJ(1));
	t108 = t106 * t96;
	t107 = t106 * t98;
	t102 = pkin(9) + qJ(3);
	t104 = cos(pkin(9));
	t103 = sin(pkin(9));
	t100 = qJ(5) + pkin(3) + pkin(8);
	t99 = cos(t102);
	t97 = sin(t102);
	t95 = sin(pkin(10)) * pkin(5) + qJ(4);
	t94 = cos(pkin(10)) * pkin(5) + qJ(2) + pkin(4) + pkin(7);
	t93 = (t100 * t104 + t103 * t95) * cos(qJ(3)) + (-t103 * t100 + t95 * t104) * sin(qJ(3)) + t104 * pkin(2) + pkin(1);
	t1 = [t97 * t108 + t109, t97 * t107 - t110, t106 * t99, t94 * t105 + t93 * t106 + 0; t97 * t110 - t107, t97 * t109 + t108, t105 * t99, t93 * t105 - t94 * t106 + 0; -t99 * t96, -t99 * t98, t97, t100 * t97 - t99 * qJ(4) + t103 * pkin(2) + 0 + pkin(6) + (sin(-pkin(10) + t102) / 0.2e1 - sin(pkin(10) + t102) / 0.2e1) * pkin(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end