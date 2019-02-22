% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP4
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
% Datum: 2019-02-22 10:40
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:40:35
% EndTime: 2019-02-22 10:40:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t81 = qJ(4) + pkin(10);
t77 = sin(t81);
t82 = sin(qJ(1));
t89 = t82 * t77;
t80 = pkin(9) + qJ(3);
t78 = cos(t80);
t88 = t82 * t78;
t79 = cos(t81);
t87 = t82 * t79;
t83 = cos(qJ(1));
t86 = t83 * t77;
t85 = t83 * t78;
t84 = t83 * t79;
t76 = sin(t80);
t75 = t78 * t84 + t89;
t74 = -t77 * t85 + t87;
t73 = -t78 * t87 + t86;
t72 = t77 * t88 + t84;
t1 = [t73, 0, -t76 * t84, t74, 0, 0; t75, 0, -t76 * t87, -t72, 0, 0; 0, 0, t78 * t79, -t76 * t77, 0, 0; t72, 0, t76 * t86, -t75, 0, 0; t74, 0, t76 * t89, t73, 0, 0; 0, 0, -t78 * t77, -t76 * t79, 0, 0; -t82 * t76, 0, t85, 0, 0, 0; t83 * t76, 0, t88, 0, 0, 0; 0, 0, t76, 0, 0, 0;];
JR_rot  = t1;
