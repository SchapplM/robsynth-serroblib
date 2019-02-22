% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:38:31
% EndTime: 2019-02-22 10:38:31
% DurationCPUTime: 0.04s
% Computational Cost: add. (62->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t90 = qJ(4) + pkin(10);
t86 = sin(t90);
t92 = sin(qJ(3));
t97 = t92 * t86;
t88 = cos(t90);
t96 = t92 * t88;
t93 = cos(qJ(3));
t95 = t93 * t86;
t94 = t93 * t88;
t91 = qJ(1) + pkin(9);
t89 = cos(t91);
t87 = sin(t91);
t85 = t87 * t86 + t89 * t94;
t84 = -t87 * t88 + t89 * t95;
t83 = -t89 * t86 + t87 * t94;
t82 = -t87 * t95 - t89 * t88;
t1 = [-t83, 0, -t89 * t96, -t84, 0, 0; t85, 0, -t87 * t96, t82, 0, 0; 0, 0, t94, -t97, 0, 0; -t87 * t92, 0, t89 * t93, 0, 0, 0; t89 * t92, 0, t87 * t93, 0, 0, 0; 0, 0, t92, 0, 0, 0; t82, 0, -t89 * t97, t85, 0, 0; t84, 0, -t87 * t97, t83, 0, 0; 0, 0, t95, t96, 0, 0;];
JR_rot  = t1;
