% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:40
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:10
	% EndTime: 2020-11-04 21:40:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:10
	% EndTime: 2020-11-04 21:40:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:10
	% EndTime: 2020-11-04 21:40:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t61 = qJ(1) + pkin(10);
	t60 = cos(t61);
	t59 = sin(t61);
	t1 = [t60, -t59, 0, cos(qJ(1)) * pkin(1) + 0; t59, t60, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:11
	% EndTime: 2020-11-04 21:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t66 = cos(qJ(3));
	t65 = sin(qJ(3));
	t64 = qJ(1) + pkin(10);
	t63 = cos(t64);
	t62 = sin(t64);
	t1 = [t63 * t66, -t63 * t65, t62, t63 * pkin(2) + t62 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t62 * t66, -t62 * t65, -t63, t62 * pkin(2) - t63 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t65, t66, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:11
	% EndTime: 2020-11-04 21:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t70 = sin(qJ(3));
	t71 = cos(qJ(3));
	t72 = pkin(3) * t71 + qJ(4) * t70 + pkin(2);
	t69 = qJ(1) + pkin(10);
	t68 = cos(t69);
	t67 = sin(t69);
	t1 = [t67, -t68 * t71, t68 * t70, cos(qJ(1)) * pkin(1) + t67 * pkin(7) + 0 + t72 * t68; -t68, -t67 * t71, t67 * t70, sin(qJ(1)) * pkin(1) - t68 * pkin(7) + 0 + t72 * t67; 0, -t70, -t71, t70 * pkin(3) - t71 * qJ(4) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:11
	% EndTime: 2020-11-04 21:40:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->21), mult. (39->24), div. (0->0), fcn. (52->8), ass. (0->13)
	t84 = pkin(3) + pkin(8);
	t83 = pkin(4) + pkin(7);
	t76 = sin(qJ(5));
	t77 = sin(qJ(3));
	t82 = t76 * t77;
	t78 = cos(qJ(5));
	t81 = t77 * t78;
	t79 = cos(qJ(3));
	t80 = qJ(4) * t77 + t84 * t79 + pkin(2);
	t75 = qJ(1) + pkin(10);
	t74 = cos(t75);
	t73 = sin(t75);
	t1 = [t73 * t78 + t74 * t82, -t73 * t76 + t74 * t81, t74 * t79, cos(qJ(1)) * pkin(1) + 0 + t83 * t73 + t80 * t74; t73 * t82 - t74 * t78, t73 * t81 + t74 * t76, t73 * t79, sin(qJ(1)) * pkin(1) + 0 - t83 * t74 + t80 * t73; -t79 * t76, -t79 * t78, t77, -t79 * qJ(4) + t84 * t77 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:11
	% EndTime: 2020-11-04 21:40:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (67->24), mult. (49->26), div. (0->0), fcn. (62->10), ass. (0->15)
	t101 = pkin(3) + pkin(9) + pkin(8);
	t100 = pkin(7) + cos(qJ(5)) * pkin(5) + pkin(4);
	t91 = qJ(5) + qJ(6);
	t88 = sin(t91);
	t93 = sin(qJ(3));
	t99 = t88 * t93;
	t89 = cos(t91);
	t98 = t89 * t93;
	t97 = pkin(5) * sin(qJ(5)) + qJ(4);
	t94 = cos(qJ(3));
	t96 = t101 * t94 + t97 * t93 + pkin(2);
	t90 = qJ(1) + pkin(10);
	t87 = cos(t90);
	t86 = sin(t90);
	t1 = [t86 * t89 + t87 * t99, -t86 * t88 + t87 * t98, t87 * t94, cos(qJ(1)) * pkin(1) + 0 + t100 * t86 + t96 * t87; t86 * t99 - t87 * t89, t86 * t98 + t87 * t88, t86 * t94, sin(qJ(1)) * pkin(1) + 0 - t100 * t87 + t96 * t86; -t94 * t88, -t94 * t89, t93, t101 * t93 - t97 * t94 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end