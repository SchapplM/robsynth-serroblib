% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:29:25
% EndTime: 2019-02-22 09:29:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (62->16), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->24)
t104 = sin(pkin(10));
t105 = sin(pkin(6));
t113 = t104 * t105;
t107 = cos(pkin(10));
t112 = t105 * t107;
t108 = cos(pkin(6));
t103 = sin(pkin(11));
t106 = cos(pkin(11));
t109 = sin(qJ(2));
t110 = cos(qJ(2));
t111 = t110 * t103 + t109 * t106;
t96 = t111 * t108;
t97 = t109 * t103 - t110 * t106;
t90 = -t104 * t97 + t107 * t96;
t92 = -t104 * t96 - t107 * t97;
t102 = qJ(4) + pkin(12);
t101 = cos(t102);
t100 = sin(t102);
t95 = t97 * t108;
t94 = t111 * t105;
t93 = t97 * t105;
t91 = t104 * t95 - t107 * t111;
t89 = -t104 * t111 - t107 * t95;
t1 = [0, t91 * t101, 0, -t92 * t100 + t101 * t113, 0, 0; 0, t89 * t101, 0, -t90 * t100 - t101 * t112, 0, 0; 0, -t93 * t101, 0, -t94 * t100 + t108 * t101, 0, 0; 0, -t91 * t100, 0, -t100 * t113 - t92 * t101, 0, 0; 0, -t89 * t100, 0, t100 * t112 - t90 * t101, 0, 0; 0, t93 * t100, 0, -t108 * t100 - t94 * t101, 0, 0; 0, t92, 0, 0, 0, 0; 0, t90, 0, 0, 0, 0; 0, t94, 0, 0, 0, 0;];
JR_rot  = t1;
