% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t63 = qJ(1) + pkin(9);
	t62 = cos(t63);
	t61 = sin(t63);
	t1 = [t62, -t61, 0, cos(qJ(1)) * pkin(1) + 0; t61, t62, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t68 = cos(qJ(3));
	t67 = sin(qJ(3));
	t66 = qJ(1) + pkin(9);
	t65 = cos(t66);
	t64 = sin(t66);
	t1 = [t65 * t68, -t65 * t67, t64, t65 * pkin(2) + t64 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t64 * t68, -t64 * t67, -t65, t64 * pkin(2) - t65 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t67, t68, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t72 = sin(qJ(3));
	t73 = cos(qJ(3));
	t74 = pkin(3) * t73 + qJ(4) * t72 + pkin(2);
	t71 = qJ(1) + pkin(9);
	t70 = cos(t71);
	t69 = sin(t71);
	t1 = [t69, -t70 * t73, t70 * t72, cos(qJ(1)) * pkin(1) + t69 * pkin(7) + 0 + t74 * t70; -t70, -t69 * t73, t69 * t72, sin(qJ(1)) * pkin(1) - t70 * pkin(7) + 0 + t74 * t69; 0, -t72, -t73, t72 * pkin(3) - t73 * qJ(4) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->21), mult. (39->24), div. (0->0), fcn. (52->8), ass. (0->13)
	t86 = pkin(3) + pkin(8);
	t85 = pkin(4) + pkin(7);
	t78 = sin(qJ(5));
	t79 = sin(qJ(3));
	t84 = t78 * t79;
	t80 = cos(qJ(5));
	t83 = t79 * t80;
	t81 = cos(qJ(3));
	t82 = qJ(4) * t79 + t86 * t81 + pkin(2);
	t77 = qJ(1) + pkin(9);
	t76 = cos(t77);
	t75 = sin(t77);
	t1 = [t75 * t80 + t76 * t84, -t75 * t78 + t76 * t83, t76 * t81, cos(qJ(1)) * pkin(1) + 0 + t85 * t75 + t82 * t76; t75 * t84 - t76 * t80, t75 * t83 + t76 * t78, t75 * t81, sin(qJ(1)) * pkin(1) + 0 - t85 * t76 + t82 * t75; -t81 * t78, -t81 * t80, t79, -t81 * qJ(4) + t86 * t79 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:36:49
	% EndTime: 2020-11-04 21:36:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (64->27), mult. (59->30), div. (0->0), fcn. (76->8), ass. (0->17)
	t102 = pkin(3) + pkin(8);
	t101 = pkin(4) + pkin(7);
	t94 = sin(qJ(5));
	t95 = sin(qJ(3));
	t100 = t94 * t95;
	t96 = cos(qJ(5));
	t99 = t95 * t96;
	t97 = cos(qJ(3));
	t98 = qJ(4) * t95 + t102 * t97 + pkin(2);
	t93 = qJ(1) + pkin(9);
	t92 = cos(t93);
	t91 = sin(t93);
	t90 = t91 * t100 - t92 * t96;
	t89 = t91 * t99 + t92 * t94;
	t88 = t92 * t100 + t91 * t96;
	t87 = t91 * t94 - t92 * t99;
	t1 = [t88, t92 * t97, t87, cos(qJ(1)) * pkin(1) + t88 * pkin(5) + t87 * qJ(6) + 0 + t101 * t91 + t98 * t92; t90, t91 * t97, -t89, sin(qJ(1)) * pkin(1) + t90 * pkin(5) - t89 * qJ(6) + 0 - t101 * t92 + t98 * t91; -t97 * t94, t95, t97 * t96, qJ(2) + pkin(6) + 0 + t102 * t95 + (-pkin(5) * t94 + qJ(6) * t96 - qJ(4)) * t97; 0, 0, 0, 1;];
	Tc_mdh = t1;
end