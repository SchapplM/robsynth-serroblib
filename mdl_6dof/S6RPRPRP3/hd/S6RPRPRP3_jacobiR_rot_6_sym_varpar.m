% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:27:03
% EndTime: 2019-02-22 10:27:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (62->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t87 = pkin(10) + qJ(5);
t83 = sin(t87);
t89 = sin(qJ(3));
t94 = t89 * t83;
t85 = cos(t87);
t93 = t89 * t85;
t90 = cos(qJ(3));
t92 = t90 * t83;
t91 = t90 * t85;
t88 = qJ(1) + pkin(9);
t86 = cos(t88);
t84 = sin(t88);
t82 = t84 * t83 + t86 * t91;
t81 = -t84 * t85 + t86 * t92;
t80 = -t86 * t83 + t84 * t91;
t79 = -t84 * t92 - t86 * t85;
t1 = [-t80, 0, -t86 * t93, 0, -t81, 0; t82, 0, -t84 * t93, 0, t79, 0; 0, 0, t91, 0, -t94, 0; -t84 * t89, 0, t86 * t90, 0, 0, 0; t86 * t89, 0, t84 * t90, 0, 0, 0; 0, 0, t89, 0, 0, 0; t79, 0, -t86 * t94, 0, t82, 0; t81, 0, -t84 * t94, 0, t80, 0; 0, 0, t92, 0, t93, 0;];
JR_rot  = t1;
